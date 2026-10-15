#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) > 0) args[which(args == "--input") + 1] else "kaks_summary.csv"

df_raw <- fread(input_file, sep = ";", fill = TRUE)

gatunki_cols <- setdiff(names(df_raw), c("gene", "parameter", "Functional_Category"))

df_wide_stats <- df_raw %>%
  select(gene, parameter, Functional_Category, all_of(gatunki_cols)) %>%
  pivot_longer(cols = all_of(gatunki_cols), names_to = "gatunek", values_to = "wartosc") %>%
  mutate(wartosc = as.numeric(wartosc)) %>%
  pivot_wider(names_from = parameter, values_from = wartosc)

if ("Ka" %in% names(df_wide_stats) && "Ks" %in% names(df_wide_stats)) {
  df_wide_stats <- df_wide_stats %>%
    drop_na(Ka, Ks) %>%
    mutate(
      Ka_Ks = Ka / Ks,
      Ka_Ks = ifelse(is.finite(Ka_Ks), Ka_Ks, NA)
    ) %>%
    drop_na(Ka_Ks)
}

max_kaks <- max(df_wide_stats$Ka_Ks, na.rm = TRUE)
if (is.infinite(max_kaks) || is.na(max_kaks) || max_kaks == 0) max_kaks <- 1.0

pA <- ggplot(df_wide_stats, aes(x = Ks, y = Ka, color = Functional_Category, size = Ka_Ks)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_size_continuous(range = c(1, 5), name = "Ka/Ks ratio") +
  theme_bw(base_size = 11) +
  labs(title = "", x = "Ks", y = "Ka", color = "Functional Category") +
  theme(legend.position = "right")

pdf("pA.pdf", width = 10, height = 6)
print(pA)
dev.off()

pB <- ggplot(df_wide_stats, aes(x = Ka_Ks, y = Functional_Category, fill = Functional_Category)) +
  geom_density_ridges(alpha = 0.7, scale = 0.9, rel_min_height = 0.01) +
  xlim(0, max_kaks * 1.05) +
  theme_bw(base_size = 11) +
  labs(title = "", x = "Ka/Ks ratio", y = "Functional Category") +
  theme(legend.position = "none")

pdf("pB.pdf", width = 10, height = 6)
print(pB)
dev.off()

pC <- ggplot(df_wide_stats, aes(x = Ka_Ks, color = Functional_Category, fill = Functional_Category)) +
  geom_density(alpha = 0.2, linewidth = 1) +
  geom_vline(xintercept = c(0.5, 1.0), linetype = "dotted", color = "red", linewidth = 0.8) +
  xlim(0, max_kaks * 1.05) +
  theme_bw(base_size = 11) +
  labs(title = "", x = "Ka/Ks ratio", y = "Density", color = "Category", fill = "Category") +
  theme(legend.position = "right")

pdf("pC.pdf", width = 10, height = 6)
print(pC)
dev.off()

selection_colors <- c(
  "Purifying (Ka/Ks < 0.5)" = "dodgerblue3", 
  "Relaxed purifying (0.5 <= Ka/Ks < 1.0)" = "goldenrod1", 
  "Neutral (Ka/Ks = 1.0)" = "gray60", 
  "Positive (Ka/Ks > 1.0)" = "firebrick2"
)

df_bar <- df_wide_stats %>%
  dplyr::mutate(
    selection_type = dplyr::case_when(
      Ka_Ks < 0.5 ~ "Purifying (Ka/Ks < 0.5)",
      Ka_Ks >= 0.5 & Ka_Ks < 1.0 ~ "Relaxed purifying (0.5 <= Ka/Ks < 1.0)",
      Ka_Ks == 1.0 ~ "Neutral (Ka/Ks = 1.0)",
      Ka_Ks > 1.0 ~ "Positive (Ka/Ks > 1.0)"
    ),
    selection_type = factor(selection_type, levels = c(
      "Purifying (Ka/Ks < 0.5)",
      "Relaxed purifying (0.5 <= Ka/Ks < 1.0)",
      "Neutral (Ka/Ks = 1.0)",
      "Positive (Ka/Ks > 1.0)"
    ))
  ) %>%
  tidyr::drop_na(selection_type) %>%
  dplyr::count(selection_type, name = "Count") %>%
  dplyr::mutate(
    prop = Count / sum(Count) * 100,
    Label = base::paste0("n = ", Count, "\n(", base::round(prop, 1), "%)")
  )

pdf("pD.pdf", width = 10, height = 6)
pD <- ggplot(df_bar, aes(x = selection_type, y = Count, fill = selection_type)) +
  geom_col(color = NA, width = 0.6) + 
  geom_text(aes(label = Label), hjust = -0.2, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = selection_colors, drop = FALSE) +
  theme_bw(base_size = 11) +
  labs(title = "", x = "Selection Category", y = "Number of Genes") +
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "bold"),
    plot.margin = margin(1, 1, 1.5, 1, "cm")
  ) +
  scale_y_continuous(
    limits = c(0, max(df_bar$Count, na.rm = TRUE) * 1.25),
    expand = c(0, 0)
  ) +
  coord_flip()
print(pD)
dev.off()
