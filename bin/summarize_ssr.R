#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

csv_files <- list.files(pattern = "_ssr_results\\.csv$")

dt_list <- list()
for (f in csv_files) {
  s_name <- sub("_ssr_results\\.csv$", "", f)
  dt_temp <- tryCatch(fread(f, sep = ";"), error = function(e) NULL)
  if (!is.null(dt_temp) && nrow(dt_temp) > 0) {
    dt_temp[, Species := s_name]
    dt_list[[f]] <- dt_temp
  }
}

full_ssr_table <- if (length(dt_list) > 0) rbindlist(dt_list, fill = TRUE) else data.table(
  ID = character(), Class = character(), Base = character(), 
  Start = numeric(), End = numeric(), Len = numeric(), 
  Count = numeric(), SSR = character(), Gene_Context = character(), 
  Species = character()
)

fwrite(full_ssr_table, "all_ssr_summary.csv", sep = ";")
cat("Aggregated SSR summary saved to all_ssr_summary.csv with", nrow(full_ssr_table), "total rows from", length(dt_list), "samples.\n")

ssr_colors <- c(
  "Mono" = "#E41A1C", 
  "Di" = "#377EB8", 
  "Tri" = "#4DAF4A", 
  "Tetra" = "#984EA3",
  "Penta" = "#FF7F00",
  "Hexa" = "#FFFF33"
)

if (nrow(full_ssr_table) > 0) {
  pdf("ssr_1.pdf", width = 6, height = 4)
  p1 <- ggplot(full_ssr_table, aes(x = Species, fill = Class)) +
    geom_bar(position = "stack") + 
    coord_flip() +                 
    scale_fill_manual(values = ssr_colors, drop = FALSE) + 
    labs(title = "", x = "Species", y = "Count", fill = "SSR Class") + 
    theme_bw() + 
    theme(axis.text.y = element_text(face = "italic"), legend.position = "right")
  print(p1)
  dev.off()
  
  pdf("ssr_2.pdf", width = 6, height = 4)
  p2 <- ggplot(full_ssr_table, aes(x = Species, y = Len, color = Class)) +
    geom_jitter(width = 0.2, size = 2.5, alpha = 0.8) +
    coord_flip() +
    scale_color_manual("SSR Class", values = ssr_colors, drop = FALSE) + 
    labs(title = "", x = "Species", y = "SSR Length (bp)") +
    theme_bw() + 
    theme(axis.text.y = element_text(face = "italic"), legend.position = "right")
  print(p2)
  dev.off()
  
  cat("SSR summary PDF plots generated successfully.\n")
} else {
  pdf("ssr_1.pdf", width = 5, height = 4)
  plot.new(); text(0.5, 0.5, "No data available")
  dev.off()
  
  pdf("ssr_2.pdf", width = 5, height = 4)
  plot.new(); text(0.5, 0.5, "No data available")
  dev.off()
}
