test_that("species spellings map to the name the reference carries", {
  # This name selects <species>/vdj/ and the <species>_<ig|tr>_<v|d|j|c> databases.
  # IgBLAST's organism vocabulary is mapped separately, where IgBLAST is called.
  expect_equal(unname(get_valid_species(words = "Homo sapiens")), "human")
  expect_equal(unname(get_valid_species(words = "MM")), "mouse")
  for (spelling in c("rhesus", "rhesus_monkey", "Rhesus macaque", "Macaca mulatta")) {
    expect_equal(unname(get_valid_species(words = spelling)), "rhesus", info = spelling)
  }
  expect_true(is.na(suppressWarnings(get_valid_species(words = "gorilla"))))
})

test_that("an unmapped species stops the input validation", {
  sheet <- tempfile(fileext = ".tsv")
  utils::write.table(data.frame(
    filename = "x.tsv", sample_id = "s1", subject_id = "p1",
    pcr_target_locus = "IGH", single_cell = "FALSE", species = "gorilla",
    tissue = "NA", sex = "NA", age = "NA", biomaterial_provider = "NA",
    stringsAsFactors = FALSE
  ), sheet, sep = "\t", row.names = FALSE, quote = FALSE)
  miairr <- system.file("rstudio/templates/project/input_validation_project_files",
                        "mapping_MiAIRR_BioSample_v1.3.1.tsv", package = "enchantr")
  skip_if(miairr == "")

  # Left to continue, the species becomes NA and reaches IgBLAST as '--organism NA'.
  expect_error(
    suppressWarnings(validate_input(sheet, miairr = miairr, collapseby = "sample_id",
                                    cloneby = "subject_id", reassign = TRUE)),
    "have no mapping"
  )
})
