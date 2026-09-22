# Build an IGH + IGL fixture for the genotype report's multi-locus path, from
# subject1 of the immcantation training data `BCR_data.tsv`, in the same columns
# as data_to_test_novel_alleles.tsv.
#
# Re-run from the package root:
#   Rscript tests/data-tests/combined_light_heavy/build_combined_light_heavy.R [path/to/BCR_data.tsv]

args <- commandArgs(trailingOnly = TRUE)
src <- if (length(args) >= 1) {
  args[[1]]
} else {
  "/home/ayelet/kleinstein_lab/immcantation/training/assets/zenodo/BCR_data.tsv"
}

out_dir <- "tests/data-tests/combined_light_heavy"
out_file <- file.path(out_dir, "combined_igh_igl.tsv.gz")

# Column contract of the heavy fixture; the combined file mirrors it exactly.
target_cols <- c(
  "sequence_id", "sequence", "rev_comp", "productive", "vj_in_frame",
  "stop_codon", "v_call", "d_call", "j_call", "sequence_alignment",
  "germline_alignment", "junction", "junction_aa", "v_cigar", "d_cigar",
  "j_cigar", "v_germline_end", "j_germline_end", "junction_length",
  "np1_length", "np2_length", "mutated_invariant", "indels", "d_5_trim",
  "d_3_trim", "j_5_trim", "subject_id", "sample_id"
)

# The training TSV quotes its string fields; let read.delim strip the quotes.
db <- read.delim(src, sep = "\t", stringsAsFactors = FALSE, quote = "\"")

# Filter on the v_call prefix; a few rows have a mismatched `locus`.
keep <- db$subject_id == "subject1" &
  !is.na(db$v_call) & nzchar(db$v_call) &
  substr(db$v_call, 1, 4) %in% c("IGHV", "IGLV")
db <- db[keep, , drop = FALSE]

# Columns the source lacks are NA.
for (col in target_cols) {
  if (!col %in% names(db)) db[[col]] <- NA
}
db$sample_id <- "sample_A"
db$sequence_id <- paste0("combined_", seq_len(nrow(db)))

db <- db[, target_cols, drop = FALSE]

# heavy then light
locus_order <- ifelse(substr(db$v_call, 1, 3) == "IGH", 1L, 2L)
db <- db[order(locus_order, db$v_call), , drop = FALSE]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
con <- gzfile(out_file, "w")
write.table(db, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
close(con)

cat("wrote", out_file, "-", nrow(db), "rows\n")
print(table(substr(db$v_call, 1, 3)))
