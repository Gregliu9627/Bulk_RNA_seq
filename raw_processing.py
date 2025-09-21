"""Raw data processing utilities for the bulk RNA-seq workflow."""
from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd

from . import config


@dataclass
class Sample:
    """Representation of a sequencing sample."""

    name: str
    fastq1: Path
    fastq2: Optional[Path] = None
    metadata: Dict[str, str] = field(default_factory=dict)

    @property
    def is_paired(self) -> bool:
        return self.fastq2 is not None


@dataclass
class RawProcessingConfig:
    """Configuration parameters for the raw processing pipeline."""

    dataset: str
    threads: int = config.DEFAULT_THREADS
    run_fastp: bool = True
    aligner: str = "salmon"
    library_type: str = "A"
    force: bool = False


class RawProcessingError(RuntimeError):
    """Raised when a step in the raw processing pipeline fails."""


class RawProcessingPipeline:
    """Pipeline orchestrating the conversion of FASTQ files to a counts matrix."""

    def __init__(self, params: RawProcessingConfig):
        self.params = params
        self.dataset_dir = config.RAW_DATA_DIR / params.dataset
        self.output_dir = config.PROCESSED_DATA_DIR / params.dataset
        self.trim_dir = self.output_dir / "trimmed"
        self.quant_dir = self.output_dir / "quant"
        self.counts_file = self.output_dir / "counts.csv"
        self.sample_sheet = self.dataset_dir / "samplesheet.csv"
        self.metadata_file = self.dataset_dir / "metadata.csv"

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def run(self) -> Path:
        self._validate_inputs()
        samples = self._load_samples()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.trim_dir.mkdir(parents=True, exist_ok=True)
        self.quant_dir.mkdir(parents=True, exist_ok=True)

        if self.params.run_fastp:
            self._run_fastp(samples)
        else:
            samples = [self._update_sample_for_existing_trim(sample) for sample in samples]

        if self.params.aligner.lower() == "salmon":
            self._run_salmon(samples)
        else:
            raise RawProcessingError(
                f"Unsupported aligner '{self.params.aligner}'. Only 'salmon' is currently supported."
            )

        self._aggregate_salmon_counts(samples)
        return self.counts_file

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _validate_inputs(self) -> None:
        if not self.dataset_dir.exists():
            raise RawProcessingError(f"Dataset directory '{self.dataset_dir}' was not found.")
        if not self.sample_sheet.exists():
            raise RawProcessingError(
                "Each dataset must contain a 'samplesheet.csv' file describing the FASTQ files."
            )
        if not config.SALMON_INDEX.exists():
            raise RawProcessingError(
                "Salmon index was not found. Please build the index in 'data/reference/salmon_index'."
            )

    def _load_samples(self) -> List[Sample]:
        sheet = pd.read_csv(self.sample_sheet)
        expected_columns = {"sample", "fastq_1"}
        if not expected_columns.issubset(sheet.columns):
            raise RawProcessingError(
                "samplesheet.csv must include at least 'sample' and 'fastq_1' columns."
            )
        samples: List[Sample] = []
        for _, row in sheet.iterrows():
            sample_name = str(row["sample"])
            fastq1 = self.dataset_dir / str(row["fastq_1"])
            fastq2_value = row.get("fastq_2")
            fastq2 = self.dataset_dir / str(fastq2_value) if pd.notna(fastq2_value) else None
            metadata = {
                col: str(row[col])
                for col in sheet.columns
                if col not in {"sample", "fastq_1", "fastq_2"}
            }
            samples.append(Sample(name=sample_name, fastq1=fastq1, fastq2=fastq2, metadata=metadata))
        return samples

    def _run_fastp(self, samples: Iterable[Sample]) -> None:
        for sample in samples:
            trimmed_r1 = self.trim_dir / f"{sample.name}_R1.fastq.gz"
            trimmed_r2 = self.trim_dir / f"{sample.name}_R2.fastq.gz" if sample.is_paired else None
            if trimmed_r1.exists() and (trimmed_r2 is None or trimmed_r2.exists()) and not self.params.force:
                sample.fastq1 = trimmed_r1
                sample.fastq2 = trimmed_r2
                continue

            cmd = [
                "fastp",
                "-i",
                str(sample.fastq1),
                "-o",
                str(trimmed_r1),
                "-w",
                str(self.params.threads),
            ]
            if sample.is_paired and sample.fastq2 is not None:
                cmd.extend(["-I", str(sample.fastq2), "-O", str(trimmed_r2)])
            self._run_command(cmd, f"fastp trimming failed for sample {sample.name}")

            sample.fastq1 = trimmed_r1
            sample.fastq2 = trimmed_r2

    def _update_sample_for_existing_trim(self, sample: Sample) -> Sample:
        trimmed_r1 = self.trim_dir / f"{sample.name}_R1.fastq.gz"
        trimmed_r2 = self.trim_dir / f"{sample.name}_R2.fastq.gz"
        if trimmed_r1.exists():
            sample.fastq1 = trimmed_r1
        if trimmed_r2.exists():
            sample.fastq2 = trimmed_r2
        return sample

    def _run_salmon(self, samples: Iterable[Sample]) -> None:
        for sample in samples:
            sample_quant_dir = self.quant_dir / sample.name
            if sample_quant_dir.exists() and not self.params.force:
                continue
            sample_quant_dir.mkdir(parents=True, exist_ok=True)

            cmd = [
                "salmon",
                "quant",
                "-i",
                str(config.SALMON_INDEX),
                "-l",
                self.params.library_type,
                "-r" if not sample.is_paired else "-1",
                str(sample.fastq1),
                "-p",
                str(self.params.threads),
                "-o",
                str(sample_quant_dir),
            ]
            if sample.is_paired and sample.fastq2 is not None:
                cmd.extend(["-2", str(sample.fastq2)])
            self._run_command(cmd, f"Salmon quantification failed for sample {sample.name}")

    def _aggregate_salmon_counts(self, samples: Iterable[Sample]) -> None:
        quant_files = {
            sample.name: self.quant_dir / sample.name / "quant.sf" for sample in samples
        }
        for sample_name, quant_file in quant_files.items():
            if not quant_file.exists():
                raise RawProcessingError(
                    f"Quantification file '{quant_file}' for sample '{sample_name}' was not found."
                )

        counts: Dict[str, pd.Series] = {}
        for sample_name, quant_file in quant_files.items():
            quant_df = pd.read_csv(quant_file, sep="\t")
            if "Name" not in quant_df or "NumReads" not in quant_df:
                raise RawProcessingError(
                    f"File '{quant_file}' does not appear to be a Salmon quantification file."
                )
            counts[sample_name] = quant_df.set_index("Name")["NumReads"]

        counts_df = pd.DataFrame(counts)
        counts_df.index.name = "transcript"
        counts_path = self.output_dir / "transcript_counts.csv"
        counts_df.to_csv(counts_path)

        # Summarise to gene-level counts when annotation is available
        if config.GTF_FILE.exists():
            gene_map = self._load_transcript_gene_map(config.GTF_FILE)
            gene_counts = (
                counts_df.assign(gene_id=lambda df: df.index.map(gene_map.get))
                .dropna(subset=["gene_id"])
                .groupby("gene_id")
                .sum()
            )
            gene_counts.index.name = "gene"
            gene_counts.to_csv(self.counts_file)
        else:
            # Fall back to transcript-level counts when annotation is missing
            counts_df.to_csv(self.counts_file)

    def _load_transcript_gene_map(self, gtf_file: Path) -> Dict[str, str]:
        mapping: Dict[str, str] = {}
        with gtf_file.open() as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                attributes = self._parse_gtf_attributes(parts[8])
                transcript_id = attributes.get("transcript_id")
                gene_id = attributes.get("gene_id")
                if transcript_id and gene_id:
                    mapping[transcript_id] = gene_id
        return mapping

    @staticmethod
    def _parse_gtf_attributes(attributes: str) -> Dict[str, str]:
        parsed: Dict[str, str] = {}
        for item in attributes.split(";"):
            item = item.strip()
            if not item:
                continue
            if " " not in item:
                continue
            key, value = item.split(" ", 1)
            parsed[key] = value.replace('"', "").strip()
        return parsed

    @staticmethod
    def _run_command(cmd: List[str], error_message: str) -> None:
        try:
            completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:  # pragma: no cover - error path
            raise RawProcessingError(
                f"{error_message}\nCommand: {' '.join(cmd)}\nStdout: {exc.stdout}\nStderr: {exc.stderr}"
            ) from exc
        else:
            log_path = config.RESULTS_DIR / "raw_processing_logs.json"
            log_path.parent.mkdir(parents=True, exist_ok=True)
            entry = {"command": cmd, "stdout": completed.stdout, "stderr": completed.stderr}
            if log_path.exists():
                data = json.loads(log_path.read_text())
                data.append(entry)
            else:
                data = [entry]
            log_path.write_text(json.dumps(data, indent=2))


def list_available_datasets() -> List[str]:
    """Return the datasets available under the raw data directory."""

    return sorted([p.name for p in config.RAW_DATA_DIR.glob("*") if p.is_dir()])


def list_processed_datasets() -> List[str]:
    """Return datasets that already have processed counts."""

    return sorted([p.name for p in config.PROCESSED_DATA_DIR.glob("*") if (p / "counts.csv").exists()])


__all__ = [
    "RawProcessingConfig",
    "RawProcessingPipeline",
    "RawProcessingError",
    "Sample",
    "list_available_datasets",
    "list_processed_datasets",
]
