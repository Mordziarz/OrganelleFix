#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Biostrings)
  library(data.table)
})


args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, args) {
  idx <- which(args == flag)

  if (length(idx) == 0L) return(NULL)
  if (length(idx) > 1L) stop("Argument provided multiple times: ", flag)
  if (idx == length(args) || startsWith(args[idx + 1L], "--")) {
    stop("Missing value for argument: ", flag)
  }

  args[idx + 1L]
}

fasta_file <- parse_arg("--fasta", args)
gff_file <- parse_arg("--gff", args)
sample_id <- parse_arg("--sample-id", args)

if (is.null(fasta_file) || !nzchar(fasta_file)) {
  stop("Error: Missing required argument --fasta.")
}
if (!file.exists(fasta_file) || dir.exists(fasta_file)) {
  stop("Error: FASTA file not found: ", fasta_file)
}
fasta_file <- normalizePath(fasta_file, mustWork = TRUE)

if (is.null(sample_id) || !nzchar(sample_id)) {
  sample_id <- sub("\\.[^.]+$", "", basename(fasta_file))
  message("Missing --sample-id, automatically set to: ", sample_id)
}
if (!is.null(gff_file) && !nzchar(gff_file)) gff_file <- NULL

if (!is.null(gff_file)) {
  if (!file.exists(gff_file) || dir.exists(gff_file)) {
    stop("Error: BED annotation file not found: ", gff_file)
  }
  gff_file <- normalizePath(gff_file, mustWork = TRUE)
}

if (!nzchar(sample_id)) stop("Error: sample-id cannot be empty.")
if (grepl("[/\\\\]", sample_id)) {
  stop("Error: sample-id must not contain directory separators.")
}

blastn_path <- Sys.which("blastn")
if (!nzchar(blastn_path)) {
  stop("Error: blastn was not found in PATH. Install BLAST+ or add blastn to PATH.")
}
message("Using BLAST executable: ", blastn_path)


deduplicate_repeat_pairs <- function(dt) {
  if (nrow(dt) == 0L) return(dt)

  required_cols <- c(
    "Sequence", "Type", "Start1", "End1", "Start2", "End2",
    "BitScore", "Length", "Identity", "Evalue"
  )

  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop(
      "Cannot deduplicate. Missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  dt <- copy(dt)

  locus1_first <- (
    dt$Start1 < dt$Start2 |
      (dt$Start1 == dt$Start2 & dt$End1 <= dt$End2)
  )

  dt[, `:=`(
    PairStart1 = ifelse(locus1_first, Start1, Start2),
    PairEnd1   = ifelse(locus1_first, End1, End2),
    PairStart2 = ifelse(locus1_first, Start2, Start1),
    PairEnd2   = ifelse(locus1_first, End2, End1)
  )]

  dt[, PairKey := paste(
    Sequence, Type, PairStart1, PairEnd1, PairStart2, PairEnd2,
    sep = "|"
  )]

  setorder(dt, PairKey, -BitScore, -Length, -Identity, Evalue)

  before <- nrow(dt)
  dt <- unique(dt, by = "PairKey")

  if (before > nrow(dt)) {
    message(
      "Deduplication: removed ",
      before - nrow(dt),
      " duplicate record(s)."
    )
  }

  dt[, c(
    "PairStart1", "PairEnd1", "PairStart2", "PairEnd2", "PairKey"
  ) := NULL]

  dt[]
}


merge_blast_hsps <- function(dt, max_gap = 50) {
  if (nrow(dt) == 0L) return(dt)

  dt <- copy(dt)

  required_cols <- c(
    "Type", "Sequence", "Start1", "End1", "Start2", "End2",
    "Length", "Identity", "Mismatches", "GapOpens",
    "BitScore", "Evalue"
  )

  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop(
      "Cannot merge HSPs. Missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  setorder(dt, Type, Sequence, Start1, Start2)

  res_list <- list()
  result_idx <- 0L

  for (tp in unique(dt$Type)) {
    sub_dt <- dt[Type == tp]
    if (nrow(sub_dt) == 0L) next

    for (seq_id in unique(sub_dt$Sequence)) {
      pair_dt <- sub_dt[Sequence == seq_id]
      setorder(pair_dt, Start1, Start2)

      merged_rows <- list()
      curr <- pair_dt[1L]

      if (nrow(pair_dt) > 1L) {
        for (i in 2:nrow(pair_dt)) {
          nxt <- pair_dt[i]

          gap1 <- nxt$Start1 - curr$End1 - 1

          if (tp %in% c("Forward", "Complement")) {
            gap2 <- nxt$Start2 - curr$End2 - 1
          } else {

            gap2 <- curr$Start2 - nxt$End2 - 1
          }

          compatible <- (
            gap1 >= 0 &&
              gap2 >= 0 &&
              gap1 <= max_gap &&
              gap2 <= max_gap &&
              abs(gap1 - gap2) <= max_gap
          )

          if (compatible) {
            old_length <- as.numeric(curr$Length)
            next_length <- as.numeric(nxt$Length)
            combined_length <- old_length + next_length

            if (combined_length > 0) {
              curr$Identity <- (
                as.numeric(curr$Identity) * old_length +
                  as.numeric(nxt$Identity) * next_length
              ) / combined_length
            }

            curr$Length <- combined_length
            curr$Mismatches <- as.numeric(curr$Mismatches) +
              as.numeric(nxt$Mismatches)
            curr$GapOpens <- as.numeric(curr$GapOpens) +
              as.numeric(nxt$GapOpens)

            curr$BitScore <- as.numeric(curr$BitScore) +
              as.numeric(nxt$BitScore)

            curr$Evalue <- min(
              as.numeric(curr$Evalue),
              as.numeric(nxt$Evalue)
            )

            curr$Start1 <- min(curr$Start1, nxt$Start1)
            curr$End1 <- max(curr$End1, nxt$End1)

            curr$Start2 <- min(curr$Start2, nxt$Start2)
            curr$End2 <- max(curr$End2, nxt$End2)

          } else {
            merged_rows[[length(merged_rows) + 1L]] <- curr
            curr <- nxt
          }
        }
      }

      merged_rows[[length(merged_rows) + 1L]] <- curr

      result_idx <- result_idx + 1L
      res_list[[result_idx]] <- rbindlist(merged_rows, fill = TRUE)
    }
  }

  if (length(res_list) == 0L) return(dt[0L])

  rbindlist(res_list, fill = TRUE)
}


find_repeats_4types_blast <- function(
    fasta_path,
    min_len = 30,
    min_identity = 90,
    evalue = 1000,
    word_size = 7,
    min_spacer = 20,
    num_threads = 4
) {
  if (!file.exists(fasta_path)) {
    stop("FASTA file not found: ", fasta_path)
  }

  raw_data <- readDNAStringSet(fasta_path)

  if (length(raw_data) == 0L) {
    warning("FASTA contains no sequences.")
    return(data.table())
  }

  all_seq_results <- list()

  outfmt <- paste(
    "6 qseqid sseqid pident length mismatch gapopen",
    "qstart qend sstart send evalue bitscore"
  )

  for (seq_idx in seq_along(raw_data)) {
    genome_seq <- raw_data[seq_idx]
    seq_id <- names(genome_seq)

    if (is.null(seq_id) || is.na(seq_id) || !nzchar(seq_id)) {
      seq_id <- paste0("seq_", seq_idx)
      names(genome_seq) <- seq_id
    }

    n <- nchar(as.character(genome_seq[[1L]]))

    message("Analyzing sequence: ", seq_id, " (", n, " bp)...")

    if (n < min_len) next

    seqs <- list(
      Forward = genome_seq,
      Reverse = reverse(genome_seq),
      Complement = complement(genome_seq),
      Palindromic = reverseComplement(genome_seq)
    )

    query_tmp <- tempfile(fileext = ".fasta")
    writeXStringSet(genome_seq, query_tmp)

    type_results <- list()

    for (tp in names(seqs)) {
      subject_seq <- seqs[[tp]]
      names(subject_seq) <- paste0(seq_id, "_", tp)

      subject_tmp <- tempfile(fileext = ".fasta")
      blast_out <- tempfile(fileext = ".tsv")

      writeXStringSet(subject_seq, subject_tmp)

      status <- tryCatch(
        system2(
          blastn_path,
          args = c(
            "-query", shQuote(normalizePath(query_tmp)),
            "-subject", shQuote(normalizePath(subject_tmp)),
            "-task", "blastn-short",
            "-word_size", as.character(word_size),
            "-evalue", as.character(evalue),
            "-strand", "plus",
            "-dust", "no",
            "-soft_masking", "false",
            "-num_threads", as.character(num_threads),
            "-outfmt", shQuote(outfmt),
            "-out", shQuote(blast_out)
          ),
          stdout = FALSE,
          stderr = FALSE
        ),
        error = function(e) {
          warning(
            "BLAST failed for sequence ", seq_id,
            ", variant ", tp, ": ", conditionMessage(e)
          )
          127L
        }
      )

      unlink(subject_tmp)

      if (status != 0L) {
        warning(
          "BLAST returned exit status ", status,
          " for ", seq_id, "/", tp
        )
        unlink(blast_out)
        next
      }

      if (!file.exists(blast_out) ||
          is.na(file.info(blast_out)$size) ||
          file.info(blast_out)$size == 0) {
        unlink(blast_out)
        next
      }

      dt <- tryCatch(
        fread(
          blast_out,
          header = FALSE,
          sep = "\t",
          showProgress = FALSE
        ),
        error = function(e) {
          warning("Could not read BLAST output: ", conditionMessage(e))
          data.table()
        }
      )

      unlink(blast_out)

      if (nrow(dt) == 0L) next

      expected_columns <- c(
        "Query", "Subject", "Identity", "Length", "Mismatches", "GapOpens",
        "qstart", "qend", "sstart", "send", "Evalue", "BitScore"
      )

      if (ncol(dt) != length(expected_columns)) {
        warning(
          "Unexpected BLAST output column count for ",
          seq_id, "/", tp
        )
        next
      }

      setnames(dt, expected_columns)

      dt <- dt[
        Length >= min_len &
          Identity >= min_identity
      ]

      if (nrow(dt) == 0L) next

      dt[, `:=`(
        Start1 = pmin(qstart, qend),
        End1 = pmax(qstart, qend)
      )]

      if (tp %in% c("Reverse", "Palindromic")) {
        dt[, `:=`(
          Start2 = n - pmax(sstart, send) + 1,
          End2 = n - pmin(sstart, send) + 1
        )]
      } else {
        dt[, `:=`(
          Start2 = pmin(sstart, send),
          End2 = pmax(sstart, send)
        )]
      }

      dt <- dt[!(Start1 == Start2 & End1 == End2)]

      if (nrow(dt) == 0L) next

      dt[, Gap := fifelse(
        End1 < Start2,
        Start2 - End1 - 1,
        fifelse(
          End2 < Start1,
          Start1 - End2 - 1,
          -1
        )
      )]

      dt <- dt[Gap >= min_spacer]

      if (nrow(dt) == 0L) next

      dt[, `:=`(
        Type = tp,
        Sequence = seq_id
      )]

      type_results[[tp]] <- dt
    }

    unlink(query_tmp)

    if (length(type_results) == 0L) next

    seq_dt <- rbindlist(type_results, fill = TRUE)

    seq_dt <- merge_blast_hsps(seq_dt, max_gap = 50)

    if (nrow(seq_dt) == 0L) next

    seq_dt[, Gap := fifelse(
      End1 < Start2,
      Start2 - End1 - 1,
      fifelse(
        End2 < Start1,
        Start1 - End2 - 1,
        -1
      )
    )]

    seq_dt <- seq_dt[Gap >= min_spacer]

    if (nrow(seq_dt) == 0L) next

    seq_dt <- deduplicate_repeat_pairs(seq_dt)

    seq_dt[, Orientation := fifelse(
      Type %in% c("Forward", "Complement"),
      "Direct",
      "Inverted"
    )]

    all_seq_results[[seq_idx]] <- seq_dt[, .(
      Sequence, Type, Orientation,
      Start1, End1, Start2, End2,
      Length, Identity, Mismatches, GapOpens,
      Gap, Evalue, BitScore
    )]
  }

  if (length(all_seq_results) == 0L) {
    return(data.table())
  }

  final_dt <- rbindlist(all_seq_results, fill = TRUE)
  setorder(final_dt, Sequence, Type, Start1, Start2)

  final_dt
}


cat(
  "Running BLAST four-transformation repeat detection for sample:",
  sample_id, "\n"
)

repeats_dt <- find_repeats_4types_blast(fasta_file)

if (nrow(repeats_dt) > 0L && !is.null(gff_file) && file.exists(gff_file)) {
  cat("Parsing BED annotations and determining gene context...\n")

  bed_dt <- fread(
    gff_file,
    header = FALSE,
    fill = TRUE,
    showProgress = FALSE
  )

  if (ncol(bed_dt) >= 4L) {
    gff_genes <- data.table(
      ID = as.character(bed_dt$V1),
      Name = as.character(bed_dt$V4),
      start = suppressWarnings(as.numeric(bed_dt$V2)),
      end = suppressWarnings(as.numeric(bed_dt$V3))
    )
    gff_genes <- gff_genes[!is.na(start) & !is.na(end)]
  } else {
    warning("BED annotation file contains fewer than four columns. Gene contexts unavailable.")
    gff_genes <- data.table(
      ID = character(), Name = character(), start = numeric(), end = numeric()
    )
  }

  get_gene_context <- function(start_pos, end_pos, gff_subset) {
    res <- character(length(start_pos))

    if (nrow(gff_subset) == 0L) {
      return(rep("Not Provided", length(start_pos)))
    }

    for (i in seq_along(start_pos)) {
      s <- min(start_pos[i], end_pos[i])
      e <- max(start_pos[i], end_pos[i])

      if (is.na(s) || is.na(e)) {
        res[i] <- NA_character_
        next
      }

      overlapped <- gff_subset[start <= e & end >= s]

      if (nrow(overlapped) > 0L) {
        res[i] <- paste(unique(overlapped$Name), collapse = ", ")
      } else {
        left_genes <- gff_subset[end < s]
        left_gen_name <- if (nrow(left_genes) > 0L) {
          left_genes[which.max(end)]$Name
        } else {
          "Start"
        }

        right_genes <- gff_subset[start > e]
        right_gen_name <- if (nrow(right_genes) > 0L) {
          right_genes[which.min(start)]$Name
        } else {
          "End"
        }

        res[i] <- paste0(
          "Intergenic (", left_gen_name, " - ", right_gen_name, ")"
        )
      }
    }

    res
  }

  for (seq_id in unique(repeats_dt$Sequence)) {
    chrom_subset <- gff_genes[ID == seq_id]

    if (nrow(chrom_subset) == 0L) chrom_subset <- gff_genes

    seq_rows <- which(repeats_dt$Sequence == seq_id)

    repeats_dt[
      seq_rows,
      Context1 := get_gene_context(Start1, End1, chrom_subset)
    ]
    repeats_dt[
      seq_rows,
      Context2 := get_gene_context(Start2, End2, chrom_subset)
    ]
  }

} else {
  if (nrow(repeats_dt) > 0L) {
    repeats_dt[, `:=`(
      Context1 = "Not Provided",
      Context2 = "Not Provided"
    )]
  } else {
    repeats_dt <- data.table(
      Sequence = character(),
      Type = character(),
      Orientation = character(),
      Start1 = numeric(),
      End1 = numeric(),
      Start2 = numeric(),
      End2 = numeric(),
      Length = numeric(),
      Identity = numeric(),
      Mismatches = numeric(),
      GapOpens = numeric(),
      Gap = numeric(),
      Evalue = numeric(),
      BitScore = numeric(),
      Context1 = character(),
      Context2 = character()
    )
  }
}


out_name <- paste0(sample_id, "_repeats_results.csv")

write.csv2(
  repeats_dt,
  out_name,
  row.names = FALSE,
  na = ""
)

cat(
  "Repeat analysis and context mapping completed successfully.\n",
  "Results saved to: ",
  normalizePath(out_name, mustWork = FALSE),
  "\n",
  sep = ""
)
