IMGT_URL <- "https://raw.githubusercontent.com/nf-core/test-datasets/airrflow/database-cache/imgtdb_base.zip"

test_that("the macaque pre-genotype filters run and report their row counts", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  input <- normalizePath(file.path("..", "data-tests", "subj_multiple_files", "db_let_12.tsv"))
  tmp_dir <- tempfile("genotype_rm_prefilters_")
  enchantr_report("genotype",
    report_params = list(
      "input" = input,
      "imgt_db" = IMGT_URL,
      "species" = "human",
      "method" = "bayesian",
      "germline_column" = "germline_alignment_d_mask",
      "unambiguous_v_only" = TRUE,
      "unmutated_v_only" = TRUE,
      "sequenced_from_v_start" = TRUE,
      "unmutated_d_only" = TRUE,
      "outdir" = tmp_dir,
      "log" = "test_rm_prefilters_command_log"
    )
  )

  steps <- read.delim(file.path(tmp_dir, "enchantr", "tables", "tab_filters.tsv"), sep = "\t")
  rows <- setNames(steps$rows_out, steps$step)
  expect_equal(
    steps$step,
    c("productive == TRUE",
      "junction_length %% 3 == 0",
      "single V assignment",
      "no mutations across the aligned V region",
      "sequence starts at the first V nucleotide",
      "no mutations in the D region (D genotype only)")
  )

  # 27 of the 1143 in-frame productive sequences carry more than one V call, and
  # 969 of the rest are mutated somewhere in V. Of the 147 left, 51 have no single
  # d_call to compare against and 2 carry a mutated D.
  expect_equal(unname(rows[3:4]), c(1116, 147))
  expect_equal(unname(rows[6]), 94)

  # Each step receives what the one before it passed on, and never invents rows.
  expect_equal(steps$rows_in[-1], steps$rows_out[-nrow(steps)])
  expect_true(all(steps$rows_out <= steps$rows_in))
})

test_that("unmutated_v_only says which germline column is empty", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  # This input carries germline_alignment_d_mask but leaves germline_alignment
  # empty, which is what happens when germlines were created with D masked.
  input <- normalizePath(file.path("..", "data-tests", "subj_multiple_files", "db_let_12.tsv"))
  tmp_dir <- tempfile("genotype_rm_germline_missing_")
  expect_error(
    enchantr_report("genotype",
      report_params = list(
        "input" = input,
        "imgt_db" = IMGT_URL,
        "species" = "human",
        "method" = "bayesian",
        "unmutated_v_only" = TRUE,
        "outdir" = tmp_dir,
        "log" = "test_rm_germline_missing_command_log"
      )
    ),
    "germline_alignment' is empty"
  )
})
