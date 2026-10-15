#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

csv_files <- list.files(pattern = "_repeats_results\\.csv$")

dt_list <- list()
for (f in csv_files) {
  s_name <- sub("_repeats_results\\.csv$", "", f)
  dt_temp <- tryCatch(fread(f, sep = ";"), error = function(e) NULL)
  if (!is.null(dt_temp) && nrow(dt_temp) > 0) {
    dt_temp[, Species := s_name]
    dt_list[[f]] <- dt_temp
  }
}

full_repeats_table <- if (length(dt_list) > 0) rbindlist(dt_list, fill = TRUE) else data.table(
  Type = character(), Start1 = numeric(), End1 = numeric(), 
  Start2 = numeric(), End2 = numeric(), Length = numeric(), 
  Identity = numeric(), Context1 = character(), Context2 = character(), 
  Species = character()
)

fwrite(full_repeats_table, "all_repeats_summary.csv", sep = ";")
cat("Aggregated summary saved to all_repeats_summary.csv with", nrow(full_repeats_table), "total rows from", length(dt_list), "samples.\n")

repeat_colors <- c(
  "Forward" = "#E41A1C", 
  "Reverse" = "#377EB8", 
  "Complement" = "#4DAF4A", 
  "Palindromic" = "#984EA3"
)

if (nrow(full_repeats_table) > 0) {
  pdf("repeat1.pdf", width = 6, height = 4)
  p1 <- ggplot(full_repeats_table, aes(x = Species, fill = Type)) +
    geom_bar(position = "stack") + 
    coord_flip() +                 
    scale_fill_manual(values = repeat_colors) + 
    labs(title = "", x = "Species", y = "Count", fill = "Repeat Type") + 
    theme_bw() + 
    theme(axis.text.y = element_text(face = "italic"), legend.position = "right")
  print(p1)
  dev.off()
  
  pdf("repeat2.pdf", width = 6, height = 4)
  p2 <- ggplot(full_repeats_table, aes(x = Species, y = Length, color = Type)) +
    geom_jitter(width = 0.2, size = 2.5, alpha = 0.8) +
    coord_flip() +
    scale_color_manual("Repeat Type", values = repeat_colors) + 
    labs(title = "", x = "Species", y = "Length (bp)") +
    theme_bw() + 
    theme(axis.text.y = element_text(face = "italic"), legend.position = "right")
  print(p2)
  dev.off()
  
  cat("Summary PDF plots generated successfully.\n")
} else {
  pdf("repeat1.pdf", width = 5, height = 4)
  plot.new(); text(0.5, 0.5, "No data available")
  dev.off()
  
  pdf("repeat2.pdf", width = 5, height = 4)
  plot.new(); text(0.5, 0.5, "No data available")
  dev.off()
}
