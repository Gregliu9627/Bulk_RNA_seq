"""Configuration for the bulk RNA-seq analysis workflow."""
from __future__ import annotations

from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent

RAW_DATA_DIR = BASE_DIR / "data" / "raw"
PROCESSED_DATA_DIR = BASE_DIR / "data" / "processed"
REFERENCE_DIR = BASE_DIR / "data" / "reference"
RESULTS_DIR = BASE_DIR / "results"

# Locations for bundled resources
SALMON_INDEX = REFERENCE_DIR / "salmon_index"
GTF_FILE = REFERENCE_DIR / "annotation.gtf"
GENE_SET_DIR = REFERENCE_DIR / "gene_sets"

SCRIPTS_DIR = BASE_DIR / "scripts"
DESEQ2_SCRIPT = SCRIPTS_DIR / "run_deseq2.R"

# Defaults used by the web interface
DEFAULT_THREADS = 8
DEFAULT_DESIGN_FORMULA = "~ condition"
DEFAULT_CONDITION_COLUMN = "condition"
DEFAULT_CONTROL_GROUP = "control"
DEFAULT_TREATMENT_GROUP = "treatment"


__all__ = [
    "BASE_DIR",
    "RAW_DATA_DIR",
    "PROCESSED_DATA_DIR",
    "REFERENCE_DIR",
    "RESULTS_DIR",
    "SALMON_INDEX",
    "GTF_FILE",
    "GENE_SET_DIR",
    "SCRIPTS_DIR",
    "DESEQ2_SCRIPT",
    "DEFAULT_THREADS",
    "DEFAULT_DESIGN_FORMULA",
    "DEFAULT_CONDITION_COLUMN",
    "DEFAULT_CONTROL_GROUP",
    "DEFAULT_TREATMENT_GROUP",
]
