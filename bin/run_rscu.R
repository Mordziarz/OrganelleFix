#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(seqinr)
library(rstatix)
library(ggpubr)
library(patchwork)
library(forcats)
library(phylogram)
library(circlize)
library(ComplexHeatmap)
library(smplot2)
library(ggtree)
library(stats)
library(stringr)
library(RSCUcaller)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide sample arguments in format ID=path/to/CDS.fasta", call. = FALSE)
}

sample_names <- c()
cds_files <- c()

for (arg in args) {
  parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    sample_names <- c(sample_names, parts[1])
    cds_files <- c(cds_files, parts[2])
  } else {
    cds_files <- c(cds_files, arg)
    sample_names <- c(sample_names, tools::file_path_sans_ext(basename(arg)))
  }
}

samples_table <- data.frame(
  sequence_path = normalizePath(cds_files),
  sample_name = sample_names,
  stringsAsFactors = FALSE
)

write.csv2(samples_table, "rscu_samples_table.csv", row.names = FALSE)

prepare_fasta(samples_table = samples_table, file_out = "combined_cds.fasta")
rscu_res_global <- get_RSCU(merged_sequences = "combined_cds.fasta")

rscu_matrix <- get_matrix(get_RSCU_out = rscu_res_global)
write.csv(rscu_matrix, "rscu_matrix_output.csv")

rscu_matrix_plastome <- rscu_matrix[grepl("plastome", rownames(rscu_matrix), ignore.case = TRUE), , drop = FALSE]
rscu_matrix_mitogenome <- rscu_matrix[grepl("mitogenome", rownames(rscu_matrix), ignore.case = TRUE), , drop = FALSE]

write.csv(rscu_matrix_plastome, "rscu_matrix_output_plastome.csv")
write.csv(rscu_matrix_mitogenome, "rscu_matrix_output_mitogenome.csv")

rscu_res_plastome <- rscu_res_global %>% filter(grepl("plastome", Species, ignore.case = TRUE))
rscu_res_mitogenome <- rscu_res_global %>% filter(grepl("mitogenome", Species, ignore.case = TRUE))

rscu_res_plastome$Species <- gsub("_plastome", "", rscu_res_plastome$Species)
rscu_res_mitogenome$Species <- gsub("_mitogenome", "", rscu_res_mitogenome$Species)

if (nrow(rscu_res_plastome) > 0) {
  pdf("rscu_heatmap_plastome.pdf", width = 14, height = 8)
  heatmap_RSCU(get_RSCU_out = rscu_res_plastome, select = "heatmap", heatmap_color = "red_blue")
  dev.off()
}

if (nrow(rscu_res_mitogenome) > 0) {
  pdf("rscu_heatmap_mitogenome.pdf", width = 14, height = 8)
  heatmap_RSCU(get_RSCU_out = rscu_res_mitogenome, select = "heatmap", heatmap_color = "red_blue")
  dev.off()
}

for (i in seq_along(sample_names)) {
  s_name <- sample_names[i]
  s_file <- cds_files[i]
  
  hist_filename <- paste0("rscu_histogram_", s_name, ".pdf")
  pdf(hist_filename, width = 8, height = 6)
  
  tryCatch({
    single_table <- data.frame(
      sequence_path = s_file,
      sample_name = s_name,
      stringsAsFactors = FALSE
    )
    single_fasta <- paste0("temp_", s_name, "_cds.fasta")
    prepare_fasta(samples_table = single_table, file_out = single_fasta)
    
    single_res <- get_RSCU(merged_sequences = single_fasta)
    
    p <- single_histogram(get_RSCU_out = single_res, title = paste("RSCU Histogram -", s_name))
    if (!is.null(p)) print(p)
    
    if (file.exists(single_fasta)) file.remove(single_fasta)
    
  }, error = function(e) {
    warning(paste("Could not generate histogram for sample:", s_name, "-", e$message))
  })
  
  dev.off()
}

cat("RSCUcaller analysis, matrix, and plots have been generated successfully with separate histograms!\n")
