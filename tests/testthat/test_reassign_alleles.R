imgt_url <- "https://raw.githubusercontent.com/nf-core/test-datasets/airrflow/database-cache/imgtdb_base.zip"

test_that("hamming distance matrix matches TIgGER mutation counts", {
  sequences <- c("acgtn.-", "ACGTAAC")
  refs <- c(ref1 = "ACGTAAA", ref2 = "ACGTTAC")

  dist_mat <- enchantr:::.compute_hamming_dist_mat(sequences, refs)
  expected <- sapply(refs, function(ref) {
    lengths(tigger::getMutatedPositions(
      sequences,
      rep(ref, length(sequences)),
      match_instead = FALSE,
      ignored_regex = "[\\.N-]"
    ))
  })
  expect_equal(unname(dist_mat), unname(expected))

  trim_data <- data.frame(d_germline_start = c(1, 2), d_germline_end = c(4, 5))
  dist_mat <- enchantr:::.compute_hamming_dist_mat(
    sequences = c("ACGT", "CGTA"),
    ref_list = refs,
    trim_seq = TRUE,
    data = trim_data,
    germline_trim_columns = c("d_germline_start", "d_germline_end"),
    indices = seq_len(nrow(trim_data))
  )
  expect_equal(dim(dist_mat), c(2L, 2L))

  expect_equal(
    enchantr:::.compute_hamming_dist_mat(c("ABCD"), c(ref = "AXYD"), ignored_regex = "Y"),
    matrix(1L, nrow = 1L, ncol = 1L)
  )
})

# TODO: find a testing set
test_that("reassign_alleles 1:1", {
  # Input in 3 files, output in 4 file.
  skip_on_cran()
  input <- normalizePath(list.files(file.path("..", "data-tests", "subj_multiple_files"), full.names = TRUE))
  tmp_dir <- file.path(tempdir(), "reassign_alleles_1_1")
  enchantr_report("reassign_alleles",
    report_params = list(
      "input" = paste(input, collapse = ","),
      "imgt_db" = imgt_url,
      "species" = "human",
      "outputby" = "subject_id",
      "outdir" = tmp_dir,
      "log" = "test_reassign_alleles_command_log"
    )
  )

  report_dir <- file.path(tmp_dir, "enchantr")
  repertoires <- list.files(file.path(report_dir, "repertoires"), full.names = TRUE)
  expect_equal(length(repertoires), 4)
  expect_equal(
    sort(basename(repertoires)),
    sort(c(
      "P05_allele-reassignment_reassign-pass.tsv",
      "Subject_0_60_allele-reassignment_reassign-pass.tsv",
      "Subject_A_allele-reassignment_reassign-pass.tsv",
      "Subject_B_allele-reassignment_reassign-pass.tsv"
    ))
  )
  db <- read_rearrangement(file.path(
    report_dir, "repertoires",
    "P05_allele-reassignment_reassign-pass.tsv"
  ))
  expect_equal(nrow(db), 1571)

  db <- read_rearrangement(file.path(
    report_dir, "repertoires",
    "Subject_0_60_allele-reassignment_reassign-pass.tsv"
  ))
  expect_equal(nrow(db), 1429)
})


test_that("reassign_alleles_after_novel", {
  # Input in 3 files, output in 4 file.
  skip_on_cran()
  input <- normalizePath(file.path("..", "data-tests", "novel_genotype", "data_to_test_novel_alleles.tsv"))
  tmp_dir <- file.path(tempdir(), "novel_inference_reassign")
  enchantr_report("novel_allele_inference",
    report_params = list(
      "input" = input,
      "imgt_db" = imgt_url,
      "species" = "human",
      "outdir" = tmp_dir,
      "pos_range" = "1:318",
      "nproc" = 1,
      "log" = "test_allele_inference_command_log"
    )
  )

  report_dir <- file.path(tmp_dir, "enchantr")
  tables <- list.files(file.path(report_dir, "tables"), full.names = T)
  novel_alleles <- read.delim(tables, sep = "\t")
  expect_equal(nrow(novel_alleles), 6)

  tmp_dir <- file.path(tempdir(), "reassign_alleles_after_novel")
  novel_db <- file.path(report_dir, "db_novel")
  enchantr_report("reassign_alleles",
    report_params = list(
      "input" = input,
      "imgt_db" = novel_db,
      "species" = "human",
      "outputby" = "subject_id",
      "outdir" = tmp_dir,
      "segments" = "v",
      "log" = "test_reassign_alleles_command_log"
    )
  )

  report_dir <- file.path(tmp_dir, "enchantr")
  repertoires <- list.files(file.path(report_dir, "repertoires"), full.names = T)
  db <- read_rearrangement(repertoires)
  expect_equal(length(grep("_", db$v_call)), 1259)
})


test_that("reassign_alleles_after_genotype_inference", {
  # Input in 3 files, output in 4 file.
  skip_on_cran()
  input <- normalizePath(file.path("..", "data-tests", "subj_multiple_files", "db_let_12.tsv"))
  tmp_dir <- file.path(tempdir(), "tigger_bayesian_genotype_reassign")
  enchantr_report("tigger_bayesian_genotype",
    report_params = list(
      "input" = input,
      "imgt_db" = imgt_url,
      "species" = "human",
      "outdir" = tmp_dir,
      "log" = "test_allele_inference_command_log"
    )
  )

  report_dir <- file.path(tmp_dir, "enchantr")
  tmp_dir <- file.path(tempdir(), "reassign_alleles_after_genotype_inference")
  genotype_db <- file.path(report_dir, "references")

  enchantr_report("reassign_alleles",
    report_params = list(
      "input" = input,
      "imgt_db" = genotype_db,
      "species" = "human",
      "outputby" = "subject_id",
      "outdir" = tmp_dir,
      "log" = "test_reassign_alleles_command_log"
    )
  )

  report_dir <- file.path(tmp_dir, "enchantr")
  repertoires <- list.files(file.path(report_dir, "repertoires"), full.names = T)
  expect_equal(length(repertoires), 2)
  db <- read_rearrangement(repertoires)
  expect_equal(nrow(db), 1144)
  expect_true(all(!is.na(db$d_call) & nzchar(db$d_call)))
  expect_true(all(!is.na(db$j_call) & nzchar(db$j_call)))
})
