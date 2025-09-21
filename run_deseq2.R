#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(readr)
  library(tibble)
  library(jsonlite)
})

option_list <- list(
  make_option(c("--counts"), type = "character", help = "Counts matrix (genes x samples)."),
  make_option(c("--metadata"), type = "character", help = "Sample metadata table."),
  make_option(c("--design"), type = "character", help = "Design formula for DESeq2."),
  make_option(c("--condition"), type = "character", help = "Metadata column for contrasts."),
  make_option(c("--control"), type = "character", help = "Control group label."),
  make_option(c("--treatment"), type = "character", help = "Treatment group label."),
  make_option(c("--output"), type = "character", help = "Path to write DE results."),
  make_option(c("--normalized"), type = "character", help = "Path to write normalized counts."),
  make_option(c("--rlog"), type = "character", help = "Path to write rlog transformed counts."),
  make_option(c("--alpha"), type = "double", default = 0.05, help = "Significance level for results()."),
  make_option(c("--lfcShrink"), type = "character", default = NULL, help = "Shrinkage method: apeglm, ashr or normal."),
  make_option(c("--diagnostics"), type = "character", default = NULL, help = "Optional JSON file for run metadata."),
  make_option(c("--format"), type = "character", default = "csv", help = "Output format: csv or tsv."),
  make_option(c("--contrast"), type = "character", default = NULL, help = "Optional custom contrast string 'factor,treatment,control'."),
  make_option(c("--blind"), type = "logical", default = FALSE, help = "Blind rlog transformation."),
  make_option(c("--minReplicatesForReplace"), type = "integer", default = 7, help = "minReplicatesForReplace for DESeq().")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

required <- c("counts", "metadata", "design", "condition", "control", "treatment", "output", "normalized", "rlog")
missing_args <- required[sapply(required, function(x) is.null(opt[[x]]))]
if (length(missing_args) > 0) {
  stop(paste("Missing required arguments:", paste(missing_args, collapse = ", ")))
}

message("Reading counts from ", opt$counts)
counts <- read_delim(opt$counts, delim = ifelse(grepl("\\.tsv$", opt$counts), "\t", ","), col_types = cols())
counts <- column_to_rownames(counts, var = colnames(counts)[1])
counts <- as.matrix(counts)

message("Reading metadata from ", opt$metadata)
metadata <- read_delim(opt$metadata, delim = ifelse(grepl("\\.tsv$", opt$metadata), "\t", ","), col_types = cols())
metadata <- column_to_rownames(metadata, var = colnames(metadata)[1])

missing_samples <- setdiff(colnames(counts), rownames(metadata))
if (length(missing_samples) > 0) {
  stop(paste("Metadata is missing samples:", paste(missing_samples, collapse = ", ")))
}
metadata <- metadata[colnames(counts), , drop = FALSE]

message("Constructing DESeqDataSet with design ", opt$design)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = as.formula(opt$design))
dds <- DESeq(dds, minReplicatesForReplace = opt$minReplicatesForReplace)

contrast <- if (!is.null(opt$contrast)) {
  strsplit(opt$contrast, ",")[[1]]
} else {
  c(opt$condition, opt$treatment, opt$control)
}
message("Using contrast: ", paste(contrast, collapse = ", "))
res <- results(dds, contrast = contrast, alpha = opt$alpha)

if (!is.null(opt$lfcShrink)) {
  res <- lfcShrink(dds, contrast = contrast, res = res, type = opt$lfcShrink)
}

res_df <- as.data.frame(res)
res_df <- rownames_to_column(res_df, var = "gene")

if (tolower(opt$format) == "tsv") {
  write_tsv(res_df, opt$output)
} else {
  write_csv(res_df, opt$output)
}

norm_counts <- counts(dds, normalized = TRUE)
norm_counts <- rownames_to_column(as.data.frame(norm_counts), var = "gene")
if (tolower(opt$format) == "tsv") {
  write_tsv(norm_counts, opt$normalized)
} else {
  write_csv(norm_counts, opt$normalized)
}

rlog_counts <- assay(rlog(dds, blind = opt$blind))
rlog_counts <- rownames_to_column(as.data.frame(rlog_counts), var = "gene")
if (tolower(opt$format) == "tsv") {
  write_tsv(rlog_counts, opt$rlog)
} else {
  write_csv(rlog_counts, opt$rlog)
}

if (!is.null(opt$diagnostics)) {
  diag <- list(
    size_factors = sizeFactors(dds),
    dispersions = dispersions(dds),
    alpha = opt$alpha,
    contrast = contrast,
    design = opt$design
  )
  write(toJSON(diag, pretty = TRUE, auto_unbox = TRUE), opt$diagnostics)
}

message("DESeq2 analysis completed.")
