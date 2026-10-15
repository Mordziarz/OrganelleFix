#!/usr/bin/env Rscript

library(metablastr)
library(circlize)
library(Biostrings)
library(dplyr)
library(rstatix)
library(ggplot2)
library(car)
library(ggpubr)
library(transgenomes)

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, args) {
  idx <- which(args == flag)
  if (length(idx) > 0) return(args[idx + 1])
  return(NULL)
}

mt_fasta <- parse_arg("--mt-fasta", args)
pt_fasta <- parse_arg("--pt-fasta", args)
mt_bed_in <- parse_arg("--mt-bed", args)
pt_bed_in <- parse_arg("--pt-bed", args)
sample_id <- parse_arg("--sample-id", args)

format_bed <- function(bed_path, out_path) {
  if (is.null(bed_path) || !file.exists(bed_path) || file.size(bed_path) == 0) {
    file.create(out_path)
    return()
  }
  dt <- read.table(bed_path, header = FALSE, stringsAsFactors = FALSE)
  
  if (ncol(dt) >= 5) {
    dt <- dt[dt[, 5] == "gene", ]
  }
  
  dt_4 <- dt[, 1:4]
  write.table(dt_4, out_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

mt_bed_clean <- "mt_clean.bed"
pt_bed_clean <- "pt_clean.bed"

format_bed(mt_bed_in, mt_bed_clean)
format_bed(pt_bed_in, pt_bed_clean)

cat("Running transfer_function for sample:", sample_id, "\n")

transfer_function_out <- transfer_function(
  fasta_mt = mt_fasta,
  fasta_pt = pt_fasta,
  bed_mt = mt_bed_clean,
  bed_pt = pt_bed_clean,
  min_length = 40,
  evalue_cut_off = 0.000001,
  min_identity = 70,
  gene_buffer = 20,
  trans_buffer = 20
)

transfers_df <- if (is.data.frame(transfer_function_out)) {
  transfer_function_out
} else if (is.list(transfer_function_out) && "transfers" %in% names(transfer_function_out)) {
  transfer_function_out$transfers
} else {
  as.data.frame(transfer_function_out)
}

write.csv(transfers_df, "transfers_table.csv", row.names = FALSE)

pdf("transfers_plot1.pdf", width = 10, height = 8)
plot_transfers(
  transfer_function_out = transfer_function_out,
  mt_sector_col = "orange",
  pt_sector_col = "darkgreen",
  mt_to_pt_col = "firebrick1",
  pt_to_mt_col = "dodgerblue1",
  unidentified_col = "grey80"
)
dev.off()

pdf("transfers_plot2.pdf", width = 10, height = 8)
plot_transfers2(
  transfer_function_out = transfer_function_out,
  normalization = TRUE,
  chromosome_gap = 0.1,
  transparency = 0.5,
  mt_sector_col = "orange",
  pt_sector_col = "darkgreen",
  mt_to_pt_col = "firebrick1",
  pt_to_mt_col = "dodgerblue1",
  unidentified_col = "grey80"
)
dev.off()

cat("Transgenomes analysis, table, and plots generated successfully!\n")
