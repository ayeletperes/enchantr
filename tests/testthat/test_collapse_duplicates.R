#### collapseDuplicates function ####

# test sequences were group by sample_id, v_call, d_call, j_call, junction_length, productive and seq_len by default
# test num_files annotation consensus_count is summed
test_that("findDuplicates", {
  # Example data.frame
  db <- data.frame(sequence_id=LETTERS[1:9],
                   sequence_alignment=c("CCCCTGGG", "CCCCTGGN","ATCGGGGA", "NAACTGGN", "NNNCTGNN", "NAACTGNG","CCCCTGGG", "CCCCTGGN","ATCGGGGA"),
                   sample_id=c('S1','S1','S1','S2','S2','S2','S2','S2','S2'),
                   v_call=c(rep("IGHV1",6),rep("IGHV2",3)),
                   d_call=rep("d_call",9),
                   j_call=rep("IGHJ2",9),
                   junction_length=rep(48,9), 
                   productive=rep(TRUE,9),
                   consensus_count=seq(1:9),
                   stringsAsFactors=FALSE)
  
  obs <- findDuplicates(db,num_fields=c('consensus_count'))
  
  exp <- data.frame(
    seq_len=rep(8,9),
    sample_id=c('S1','S1','S1','S2','S2','S2','S2','S2','S2'),
    sequence_id=LETTERS[1:9],
    sequence_alignment=c("CCCCTGGG", "CCCCTGGN","ATCGGGGA", "NAACTGGN", "NNNCTGNN", "NAACTGNG","CCCCTGGG", "CCCCTGGN","ATCGGGGA"),
    v_call=c(rep("IGHV1",6),rep("IGHV2",3)),
    d_call=rep("d_call",9),
    j_call=rep("IGHJ2",9),
    junction_length=rep(48,9),
    productive=rep(TRUE,9),
    collapse_pass=c(TRUE,FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE),
    collapse_count=c(2,NA,1,3,NA,NA,2,NA,1),
    consensus_count=c(3,NA,3,15,NA,NA,15,NA,9),
    stringsAsFactors = F
  )
  
  expect_equivalent(obs, exp)
})

test_that("findDuplicates respects cell_id in single-cell data", {
  # Two records from different cells share the same sequence — they must NOT be collapsed.
  # One additional record from cell2 has a unique sequence and should pass as-is.
  db <- data.frame(
    sequence_id        = c("cell1_contig1", "cell2_contig1", "cell2_contig2"),
    sequence_alignment = c("CCCCTGGG",      "CCCCTGGG",      "ATCGGGGA"),
    cell_id            = c("cell1",         "cell2",         "cell2"),
    sample_id          = rep("S1", 3),
    v_call             = rep("IGHV1", 3),
    d_call             = rep("d_call", 3),
    j_call             = rep("IGHJ2", 3),
    junction_length    = rep(48, 3),
    productive         = rep(TRUE, 3),
    stringsAsFactors   = FALSE
  )

  obs <- findDuplicates(db, groups = "sample_id", num_fields = NULL)

  # All three records should pass: identical sequences from different cells must be kept separate.
  expect_true(all(obs$collapse_pass))
  expect_equal(nrow(obs[obs$collapse_pass, ]), 3)
})



test_that("mask_5prime_sequence_alignment", {
  # Example data.frame
  db <- data.frame(sequence_id=c('A'),
                   sequence_alignment=c(".................................................CCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGAGGCACCTTC............AGCAGCTATGCTATCAGCTGGGTGCGACAGGCCCCCGGACAAGGGCTTGAGTGGATGGGAAGGATCATCCCTATC......CTTGGTATAGCAAACTACGCACAGAAGTTCCAG...GGCAGAGTCACGATTACCGCGGACAAATCCACGAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTNTGCGAGTCACTACGATTTTTGGAGTGGTTATCCCCCGGGATACTACTACTACGGTATGGACGTCTGGGGCCAAGGAACCACGGTCACCGTCTCCTCAG"),
                   v_germline_end=c(317),
                   stringsAsFactors=FALSE)
  
  obs <- mask_5prime_sequence_alignment(db,80)
  
  exp <- data.frame(
    sequence_id=c('A'),
    sequence_alignment=c(".................................................NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGGCACCTTC............AGCAGCTATGCTATCAGCTGGGTGCGACAGGCCCCCGGACAAGGGCTTGAGTGGATGGGAAGGATCATCCCTATC......CTTGGTATAGCAAACTACGCACAGAAGTTCCAG...GGCAGAGTCACGATTACCGCGGACAAATCCACGAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTNTGCGAGTCACTACGATTTTTGGAGTGGTTATCCCCCGGGATACTACTACTACGGTATGGACGTCTGGGGCCAAGGAACCACGGTCACCGTCTCCTCAG"),
    v_germline_end=c(317),
    stringsAsFactors=FALSE
  )
  
  expect_equivalent(obs, exp)
})


test_that("mask_3prime_sequence_alignment", {
  # Example data.frame
  db <- data.frame(sequence_id=c('A','B'),
                   sequence_alignment=c(".................................................CCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGAGGCACCTTC............AGCAGCTATGCTATCAGCTGGGTGCGACAGGCCCCCGGACAAGGGCTTGAGTGGATGGGAAGGATCATCCCTATC......CTTGGTATAGCAAACTACGCACAGAAGTTCCAG...GGCAGAGTCACGATTACCGCGGACAAATCCACGAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTNTGCGAGTCACTACGATTTTTGGAGTGGTTATCCCCCGGGATACTACTACTACGGTATGGACGTCTGGGGCCAAGGAACCACGGTCACCGTCTCCTCAG", "CCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGAGGCACCTT"),
                   stringsAsFactors=FALSE)
  
  obs <- mask_3prime_sequence_alignment(db,25)
  
  exp <- data.frame(
    sequence_id=c('A','B'),
    sequence_alignment=c(".................................................CCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGAGGCACCTTC............AGCAGCTATGCTATCAGCTGGGTGCGACAGGCCCCCGGACAAGGGCTTGAGTGGATGGGAAGGATCATCCCTATC......CTTGGTATAGCAAACTACGCACAGAAGTTCCAG...GGCAGAGTCACGATTACCGCGGACAAATCCACGAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTNTGCGAGTCACTACGATTTTTGGAGTGGTTATCCCCCGGGATACTACTACTACGGTATGGACGTCTGGGGCCAANNNNNNNNNNNNNNNNNNNNNNNNN","CCTCAGTGAAGGTCTNNNNNNNNNNNNNNNNNNNNNNNNN"),
    stringsAsFactors=FALSE
  )
  
  expect_equivalent(obs, exp)
})


#### Collapse Duplicates Report ####
test_that("Collapse duplicates on sample_id", {
  # Input in one file, output in 4 files.
  skip_on_cran()
  input <- normalizePath(file.path("..", "data-tests", "subj_multiple_files","data_to_test_collapse_duplicates.tsv"))
  tmp_dir <- file.path(tempdir(),"collapse_duplicates_on_sample")
  report_params <- list(
    input = input,
    collapseby = "sample_id", outputby = "sample_id",
    outdir = tmp_dir,
    nproc = 1,
    log = "test_collapse_duplicate_command_log"
  )
  suppressWarnings(enchantr_report("collapse_duplicates", report_params = report_params))
  report_dir <- file.path(tmp_dir, "enchantr")
  repertoires <- list.files(file.path(report_dir, "repertoires"), full.names = TRUE)
  # test number and name of output files, test number of sequence and collapse_count in file S2_collapse_collapse-pass.tsv.gz
  expect_equal(length(repertoires), 4)
  expect_equal(basename(repertoires),c('S1_collapse_collapse-pass.tsv.gz','S2_collapse_collapse-pass.tsv.gz','S3_collapse_collapse-pass.tsv.gz','S4_collapse_collapse-pass.tsv.gz'))
  db_s2 <- suppressWarnings(read_rearrangement(file.path(report_dir, "repertoires", "S2_collapse_collapse-pass.tsv.gz")))
  expect_equal(nrow(db_s2), 3)
  expect_equal(db_s2$collapse_count, c('3','2','1'))
})


test_that("Collapse duplicates on subject_id", {
  # Input in one file, output in 2 files.
  skip_on_cran()
  input <- normalizePath(file.path("..", "data-tests", "subj_multiple_files","data_to_test_collapse_duplicates.tsv"))
  tmp_dir_subject <- file.path(tempdir(),"collapse_duplicates_on_subject")
  report_params <- list(
    input = input,
    collapseby = "subject_id", outputby = "subject_id",
    outdir = tmp_dir_subject,
    nproc = 1,
    log = "test_collapse_duplicate_command_log"
  )
  suppressWarnings(enchantr_report("collapse_duplicates", report_params = report_params))
  report_dir <- file.path(tmp_dir_subject, "enchantr")
  repertoires <- list.files(file.path(report_dir, "repertoires"), full.names = TRUE)
  # test number of output files, name of output files, number of sequence , collapse_count and consensus count in file Subject_B_collapse_collapse-pass.tsv
  expect_equal(length(repertoires), 2)
  expect_equal(basename(repertoires),c('Subject_A_collapse_collapse-pass.tsv.gz','Subject_B_collapse_collapse-pass.tsv.gz'))
  db_B <- suppressWarnings(read_rearrangement(file.path(report_dir, "repertoires", "Subject_B_collapse_collapse-pass.tsv.gz")))
  expect_equal(nrow(db_B), 5)
  expect_equal(db_B$collapse_count, c('1','1','1','2','4'))
  expect_equal(db_B$consensus_count, c(10,11,12,29,64))
})


test_that("Collapse duplicates on sample_id, filter out collapse_count<2", {
  # Input in one file, output in 4 files.
  skip_on_cran()
  input <- normalizePath(file.path("..", "data-tests", "subj_multiple_files","data_to_test_collapse_duplicates.tsv"))
  tmp_dir <- file.path(tempdir(),"collapse_duplicates_on_sample_filter_count")
  report_params <- list(
    'input' = input,
    'collapseby' = "sample_id", 'outputby' = "sample_id",
    'collapse_filter_threshold' = 2,
    'outdir' = tmp_dir,
    'nproc' = 1,
    'log' = "test_collapse_duplicate_command_log"
  )
  suppressWarnings(enchantr_report("collapse_duplicates", report_params = report_params))
  report_dir <- file.path(tmp_dir, "enchantr")
  repertoires <- list.files(file.path(report_dir, "repertoires"), full.names = TRUE)
  # test number and name of output files, test number of sequence and collapse_count in file S2_collapse_collapse-pass.tsv.gz
  expect_equal(length(repertoires), 4)
  expect_equal(basename(repertoires),c('S1_collapse_collapse-pass.tsv.gz','S2_collapse_collapse-pass.tsv.gz','S3_collapse_collapse-pass.tsv.gz','S4_collapse_collapse-pass.tsv.gz'))
  db_s2 <- suppressWarnings(read_rearrangement(file.path(report_dir, "repertoires", "S2_collapse_collapse-pass.tsv.gz")))
  expect_equal(nrow(db_s2), 2)
  expect_equal(db_s2$collapse_count, c('3','2'))
})


test_that("Collapse duplicates on sample_id, mask at IMGT position 3", {
  # Input in one file, output in 4 files.
  skip_on_cran()
  input <- normalizePath(file.path("..", "data-tests", "subj_multiple_files","data_to_test_collapse_duplicates.tsv"))
  tmp_dir <- file.path(tempdir(),"collapse_duplicates_on_sample_IMGT_masking")
  report_params <- list(
    'input' = input,
    'collapseby' = "sample_id", 'outputby' = "sample_id",
    'mask_imgt_position' = 3,
    'outdir' = tmp_dir,
    'nproc' = 1,
    'log' = "test_collapse_duplicate_command_log"
  )
  suppressWarnings(enchantr_report("collapse_duplicates", report_params = report_params))
  report_dir <- file.path(tmp_dir, "enchantr")
  repertoires <- list.files(file.path(report_dir, "repertoires"), full.names = TRUE)
  # test number and name of output files, test number of sequence, collapse_count and consunsus count in S3_collapse_collapse-pass.tsv.gz
  expect_equal(length(repertoires), 4)
  expect_equal(basename(repertoires),c('S1_collapse_collapse-pass.tsv.gz','S2_collapse_collapse-pass.tsv.gz','S3_collapse_collapse-pass.tsv.gz','S4_collapse_collapse-pass.tsv.gz'))
  db_s3 <- suppressWarnings(read_rearrangement(file.path(report_dir, "repertoires", "S3_collapse_collapse-pass.tsv.gz")))
  expect_equal(nrow(db_s3), 2)
  expect_equal(db_s3$collapse_count, c('4','2'))
  expect_equal(db_s3$consensus_count, c(50,25))
})


test_that("Collapse duplicates on sample_id, mask 3 bases to 3 prime end", {
  # Input in one file, output in 4 files.
  skip_on_cran()
  input <- normalizePath(file.path("..", "data-tests", "subj_multiple_files","data_to_test_collapse_duplicates.tsv"))
  tmp_dir <- file.path(tempdir(),"collapse_duplicates_on_sample_IMGT_masking")
  report_params <- list(
    'input' = input,
    'collapseby' = "sample_id", 'outputby' = "sample_id",
    'mask_length_to_3end' = 3,
    'outdir' = tmp_dir,
    'nproc' = 1,
    'log' = "test_collapse_duplicate_command_log"
  )
  suppressWarnings(enchantr_report("collapse_duplicates", report_params = report_params))
  report_dir <- file.path(tmp_dir, "enchantr")
  repertoires <- list.files(file.path(report_dir, "repertoires"), full.names = TRUE)
  # test number and name of output files, test sequence_alignment in S3_collapse_collapse-pass.tsv.gz
  expect_equal(length(repertoires), 4)
  expect_equal(basename(repertoires),c('S1_collapse_collapse-pass.tsv.gz','S2_collapse_collapse-pass.tsv.gz','S3_collapse_collapse-pass.tsv.gz','S4_collapse_collapse-pass.tsv.gz'))
  db_s3 <- suppressWarnings(read_rearrangement(file.path(report_dir, "repertoires", "S3_collapse_collapse-pass.tsv.gz")))
  expect_equal(db_s3[['sequence_alignment']], c('CCCCTNNN', 'ACCCTNNN', 'ATCGGNNN', 'CTCGGNNN','NAACTNNN'))
})

test_that("sequences that do not start at the first V nucleotide are counted, and dropped on request", {
  skip_on_cran()
  # 21 of the 1144 sequences in this repertoire align from a later position.
  input <- normalizePath(file.path("..", "data-tests", "subj_multiple_files", "db_let_12.tsv"))
  tmp_dir <- file.path(tempdir(), "collapse_duplicates_v_start")
  suppressWarnings(enchantr_report("collapse_duplicates", report_params = list(
    input = input,
    collapseby = "subject_id", outputby = "subject_id",
    sequenced_from_v_start = TRUE,
    outdir = tmp_dir,
    nproc = 1,
    log = "test_collapse_duplicates_v_start_command_log"
  )))

  report_dir <- file.path(tmp_dir, "enchantr")
  counts <- read.delim(file.path(report_dir, "tables", "tab_v_start.tsv"), sep = "\t")
  expect_equal(sum(counts$sequences), 1144)
  expect_equal(sum(counts$not_from_v_start), 21)

  # None of them survive into the output that clonal analysis reads.
  kept <- do.call(rbind, lapply(
    list.files(file.path(report_dir, "repertoires"), full.names = TRUE),
    function(f) suppressWarnings(read_rearrangement(f))
  ))
  expect_false(any(startsWith(kept$sequence_alignment, ".")))
})
