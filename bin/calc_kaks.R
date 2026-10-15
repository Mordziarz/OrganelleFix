#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(seqinr)
  library(data.table)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx + 1 <= length(args)) {
    return(args[idx + 1])
  }
  return(NULL)
}

organelle <- parse_arg("--organelle")
samples_arg <- parse_arg("--samples")
cds_arg <- parse_arg("--cds_files")
reference_arg <- parse_arg("--reference")

if (is.null(organelle) || is.null(samples_arg) || is.null(cds_arg) || is.null(reference_arg)) {
  stop("Missing required arguments: --organelle, --samples, --cds_files, or --reference")
}

samples <- unlist(strsplit(samples_arg, ","))
cds_files <- unlist(strsplit(cds_arg, ","))

if (length(samples) != length(cds_files)) {
  stop("Number of samples does not match number of CDS files!")
}

ref_sample <- reference_arg
if (!(ref_sample %in% samples)) {
  stop(sprintf("Specified reference sample '%s' not found in samples list!", ref_sample))
}

cat(sprintf("[INFO] Processing Ka/Ks calculation for organelle: %s with %d samples. Reference: %s\n", 
            organelle, length(samples), ref_sample))

get_functional_category <- function(gene_name) {
  gene <- tolower(gene_name)
  if (grepl("psa", gene)) {
    return("Photosystem I")
  } else if (grepl("psb", gene)) {
    return("Photosystem II")
  } else if (grepl("pet", gene)) {
    return("Cytochrome b/f complex")
  } else if (grepl("rbcl", gene)) {
    return("RUBISCO large subunit")
  } else if (grepl("atp", gene)) {
    return("ATP synthase")
  } else if (grepl("rpo", gene)) {
    return("Dna dependent RNA polymerase")
  } else if (grepl("rps", gene)) {
    return("Small subunit of ribosome")
  } else if (grepl("rpl", gene)) {
    return("Large subunit of ribosome")
  } else if (grepl("ndh|nad", gene)) {
    return("NADH dehydrogenase")
  } else if (grepl("cox|cob|ccm|sdh", gene)) {
    return("Mitochondrial Complex / Respiration")
  } else {
    return("Other genes")
  }
}
normalize_name <- function(header) {
  header_clean <- gsub("^>", "", header)
  
  if (grepl("gene=", header_clean, ignore.case = TRUE)) {
    m <- sub(".*gene=([a-zA-Z0-9_-]+).*", "\\1", header_clean, ignore.case = TRUE)
    return(tolower(gsub("[^a-zA-Z0-9]", "", m)))
  }
  
  parts <- unlist(strsplit(header_clean, "[_\\s\\[\\]()]+"))
  parts <- parts[parts != ""]
  
  known_patterns <- c("nad", "cox", "cob", "atp", "psa", "psb", "pet", "rpo", 
                      "rps", "rpl", "rrn", "trn", "ccm", "sdh", "matr", "ycf", 
                      "rbcl", "accd", "cema", "ccsa", "ndh")
  
  for (p in parts) {
    p_lower <- tolower(p)
    for (kp in known_patterns) {
      if (grepl(kp, p_lower)) {
        return(gsub("[^a-zA-Z0-9]", "", p_lower))
      }
    }
  }
  
  if (length(parts) > 0) {
    last_part <- parts[length(parts)]
    return(tolower(gsub("[^a-zA-Z0-9]", "", last_part)))
  }
  
  return("unnamed")
}

sample_data <- list()
all_genes <- c()

for (i in seq_along(samples)) {
  s_name <- samples[i]
  c_path <- cds_files[i]
  
  if (!file.exists(c_path)) {
    warning(paste("File not found:", c_path))
    next
  }
  
  fa <- read.fasta(c_path, seqtype = "DNA", set.attributes = TRUE)
  genes_list <- list()
  
  for (seq_name in names(fa)) {
    annot <- getAnnot(fa[[seq_name]])
    g_name <- normalize_name(annot)
    if (g_name != "" && g_name != "unnamed") {
      seq_seq <- toupper(c2s(fa[[seq_name]]))
      genes_list[[g_name]] <- seq_seq
    }
  }
  sample_data[[s_name]] <- genes_list
  all_genes <- union(all_genes, names(genes_list))
}

all_genes <- sort(all_genes)
output_rows <- list()

cat(sprintf("[INFO] Unique genes found across samples: %d\n", length(all_genes)))

for (gene in all_genes) {
  if (nchar(gene) == 0) next
  
  category <- get_functional_category(gene)
  
  if (is.null(sample_data[[ref_sample]][[gene]])) {
    next
  }
  
  ref_seq_str <- sample_data[[ref_sample]][[gene]]
  
  if (nchar(ref_seq_str) %% 3 != 0) {
    next
  }
  
  ka_row <- list(gene = gene, parameter = "Ka", Functional_Category = category)
  ks_row <- list(gene = gene, parameter = "Ks", Functional_Category = category)
  kaks_row <- list(gene = gene, parameter = "Ka/Ks", Functional_Category = category)
  
  for (s in samples) {
    if (s == ref_sample) {
      ka_row[[s]] <- 0.0
      ks_row[[s]] <- 0.0
      kaks_row[[s]] <- 0.0
      next
    }
    
    if (is.null(sample_data[[s]][[gene]])) {
      ka_row[[s]] <- NA_real_
      ks_row[[s]] <- NA_real_
      kaks_row[[s]] <- NA_real_
      next
    }
    
    s_seq_str <- sample_data[[s]][[gene]]
    
    if (nchar(s_seq_str) %% 3 != 0) {
      ka_row[[s]] <- NA_real_
      ks_row[[s]] <- NA_real_
      kaks_row[[s]] <- NA_real_
      next
    }
    
    tmp_fa <- tempfile(fileext = ".fasta")
    writeLines(c(paste0(">", ref_sample), ref_seq_str, paste0(">", s), s_seq_str), tmp_fa)
    
    tmp_aln <- paste0(tmp_fa, ".aln")
    mafft_cmd <- sprintf("mafft --auto --quiet %s > %s 2>/dev/null", tmp_fa, tmp_aln)
    exit_code <- system(mafft_cmd)
    
    if (exit_code != 0 || !file.exists(tmp_aln) || file.info(tmp_aln)$size == 0) {
      if (file.exists(tmp_fa)) file.remove(tmp_fa)
      if (file.exists(tmp_aln)) file.remove(tmp_aln)
      ka_row[[s]] <- NA_real_; ks_row[[s]] <- NA_real_; kaks_row[[s]] <- NA_real_
      next
    }
    
    aln_obj <- tryCatch({
      read.alignment(file = tmp_aln, format = "fasta")
    }, error = function(e) { NULL })
    
    if (file.exists(tmp_fa)) file.remove(tmp_fa)
    if (file.exists(tmp_aln)) file.remove(tmp_aln)
    
    if (is.null(aln_obj) || length(aln_obj$nam) < 2) {
      ka_row[[s]] <- NA_real_; ks_row[[s]] <- NA_real_; kaks_row[[s]] <- NA_real_
      next
    }
    
    # Bezpieczne wywołanie kaks() z obsługą błędów i ostrzeżeń
    kaks_res <- suppressWarnings(tryCatch({
      kaks(aln_obj)
    }, error = function(e) { NULL }))
    
    # Rygorystyczne sprawdzenie, czy wynik jest poprawną listą z macierzami ka i ks
    if (is.null(kaks_res) || !is.list(kaks_res) || is.null(kaks_res$ka) || is.null(kaks_res$ks)) {
      ka_row[[s]] <- NA_real_
      ks_row[[s]] <- NA_real_
      kaks_row[[s]] <- NA_real_
      next
    }
    
    ka_mat <- kaks_res$ka
    ks_mat <- kaks_res$ks
    
    ref_idx <- which(aln_obj$nam == ref_sample)
    s_idx <- which(aln_obj$nam == s)
    
    ka_val <- NA
    ks_val <- NA
    
    if (is.matrix(ka_mat)) {
      if (length(ref_idx) > 0 && length(s_idx) > 0 &&
          ref_idx <= nrow(ka_mat) && s_idx <= ncol(ka_mat)) {
        ka_val <- ka_mat[s_idx, ref_idx]
        ks_val <- ks_mat[s_idx, ref_idx]
      }
    } else if (length(ka_mat) == 1) {
      ka_val <- ka_mat[1]
      ks_val <- ks_mat[1]
    }
    
    if (is.na(ka_val) || is.na(ks_val) || ka_val < 0 || ks_val < 0) {
      ka_row[[s]] <- NA_real_
      ks_row[[s]] <- NA_real_
      kaks_row[[s]] <- NA_real_
    } else {
      k_a <- as.numeric(ka_val)
      k_s <- as.numeric(ks_val)
      ka_row[[s]] <- round(k_a, 5)
      ks_row[[s]] <- round(k_s, 5)
      kaks_row[[s]] <- ifelse(k_s > 0, round(k_a / k_s, 4), 0.0)
    }
  }
  
  output_rows[[length(output_rows) + 1]] <- as.data.frame(ka_row, stringsAsFactors = FALSE)
  output_rows[[length(output_rows) + 1]] <- as.data.frame(ks_row, stringsAsFactors = FALSE)
  output_rows[[length(output_rows) + 1]] <- as.data.frame(kaks_row, stringsAsFactors = FALSE)
}

if (length(output_rows) == 0) {
  warning("No Ka/Ks results could be computed!")
  empty_df <- data.frame(gene = character(), parameter = character(), Functional_Category = character())
  for (s in samples) empty_df[[s]] <- numeric()
  fwrite(empty_df, "kaks_summary.csv", sep = ";", na = "NA")
} else {
  final_dt <- bind_rows(output_rows)
  fwrite(final_dt, "kaks_summary.csv", sep = ";", na = "NA")
  cat(sprintf("[INFO] Generated Ka/Ks summary table for %s: kaks_summary.csv\n", organelle))
}
