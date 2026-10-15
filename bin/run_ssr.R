#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, args) {
  idx <- which(args == flag)
  if (length(idx) > 0) return(args[idx + 1])
  return(NULL)
}

fasta_file <- parse_arg("--fasta", args)
gff_file <- parse_arg("--gff", args)
sample_id <- parse_arg("--sample-id", args)

analyze_ssr_final <- function(fasta_path, 
                            min_reps = c(Mono=12, Di=6, Tri=4, Tetra=3, Penta=3, Hexa=3)) {
  
  lines <- readLines(fasta_path)
  header_indices <- grep("^>", lines)
  if(length(header_indices) == 0) return(data.table())
  
  all_results <- list()
  
  for (h in 1:length(header_indices)) {
    start_idx <- header_indices[h]
    end_idx <- ifelse(h < length(header_indices), header_indices[h+1] - 1, length(lines))
    
    raw_id <- sub("^>", "", lines[start_idx])
    id <- strsplit(raw_id, "\\s+")[[1]][1]
    
    full_seq <- paste(lines[(start_idx + 1):end_idx], collapse = "")
    n <- nchar(full_seq)
    
    for (m_len in 1:6) {
      class_name <- c("Mono", "Di", "Tri", "Tetra", "Penta", "Hexa")[m_len]
      req_reps <- min_reps[class_name]
      
      max_start <- n - (req_reps * m_len) + 1
      if (max_start < 1) next
      
      for (i in 1:max_start) {
        motif <- substr(full_seq, i, i + m_len - 1)
        
        count <- 1
        while(TRUE) {
          next_start <- i + (count * m_len)
          next_end <- next_start + m_len - 1
          if (next_end > n) break
          if (substr(full_seq, next_start, next_end) == motif) {
            count <- count + 1
          } else {
            break
          }
        }
        
        if (count >= req_reps) {
          all_results[[length(all_results) + 1]] <- data.table(
            ID = id,
            Class = class_name,
            Base = motif,
            Start = i,
            End = i + (count * m_len) - 1,
            Len = count * m_len,
            Count = count,
            SSR = substr(full_seq, i, i + (count * m_len) - 1)
          )
        }
      }
    }
  }
  
  if (length(all_results) == 0) return(data.table())
  
  dt <- rbindlist(all_results)
  dt <- dt[order(ID, -Len)]
  
  final_dt <- data.table()
  
  for (current_id in unique(dt$ID)) {
    subset_dt <- dt[ID == current_id]
    kept <- data.table()
    
    while(nrow(subset_dt) > 0) {
      top <- subset_dt[1]
      kept <- rbind(kept, top)
      subset_dt <- subset_dt[!(Start <= top$End & End >= top$Start)]
    }
    final_dt <- rbind(final_dt, kept)
  }
  
  return(final_dt[order(ID, Start)])
}

cat("Analyzing SSRs for sample:", sample_id, "\n")
ssr_results <- analyze_ssr_final(fasta_file)

if (nrow(ssr_results) > 0 && !is.null(gff_file) && file.exists(gff_file)) {
  cat("Parsing annotations and determining gene context...\n")
  
  bed_dt <- fread(gff_file, header = FALSE, fill = TRUE)
  
  if (ncol(bed_dt) >= 4) {
    gff3_genes <- data.table(
      ID = as.character(bed_dt$V1), 
      Name = as.character(bed_dt$V4), 
      start = as.numeric(bed_dt$V2), 
      end = as.numeric(bed_dt$V3)
    )
  } else {
    gff3_genes <- data.table(ID = character(0), Name = character(0), start = numeric(0), end = numeric(0))
  }
  
  get_gene_context <- function(start_pos, end_pos, seq_id, gff_subset) {
    res <- character(length(start_pos))
    
    for (i in seq_along(start_pos)) {
      s <- min(start_pos[i], end_pos[i])
      e <- max(start_pos[i], end_pos[i])
      curr_id <- as.character(seq_id[i])
      
      if (is.na(s) || is.na(e)) {
        res[i] <- NA
        next
      }
      
      chr_genes <- gff_subset[ID == curr_id]
      if (nrow(chr_genes) == 0) {
        res[i] <- "Intergenic (Unknown - Unknown)"
        next
      }
      
      overlapped <- chr_genes[start <= e & end >= s]
      
      if (nrow(overlapped) > 0) {
        res[i] <- paste(unique(overlapped$Name), collapse = ", ")
      } else {
        left_genes <- chr_genes[end < s]
        if (nrow(left_genes) > 0) {
          closest_left <- left_genes[which.max(end)]
          left_gen_name <- closest_left$Name
        } else {
          left_gen_name <- "Start"
        }
        
        right_genes <- chr_genes[start > e]
        if (nrow(right_genes) > 0) {
          closest_right <- right_genes[which.min(start)]
          right_gen_name <- closest_right$Name
        } else {
          right_gen_name <- "End"
        }
        
        res[i] <- paste0("Intergenic (", left_gen_name, " - ", right_gen_name, ")")
      }
    }
    return(res)
  }
  
  ssr_results[, Gene_Context := get_gene_context(Start, End, ID, gff3_genes)]
} else {
  ssr_results[, Gene_Context := "Not Provided"]
}

out_name <- paste0(sample_id, "_ssr_results.csv")
write.csv2(ssr_results, out_name, row.names = FALSE)
cat("SSR analysis and context mapping completed successfully. Results saved to:", out_name, "\n")
