db <- tibble::tibble(
  accession = c("A", "B", "C", "D", "E"),
  name = c(
    "genome_a",
    "genome_b",
    "genome_c",
    "genome_d",
    "genome_e"
  ),
  taxid = c(1, 2, 3, 4, 5),
  filename = c(
    "a.fna",
    "b.fna",
    "c.fna",
    "d.fna",
    "e.fna"
  ),
  assembly_level = c(
    "Complete Genome",
    "Chromosome",
    "Scaffold",
    "Contig",
    NA
  ),
  contig_n50 = c(
    100000,
    50000,
    20000,
    5000,
    NA
  ),
  scaffold_n50 = c(
    200000,
    100000,
    50000,
    10000,
    NA
  ),
  check_m_completeness = c(
    99,
    95,
    90,
    70,
    NA
  ),
  check_m_contamination = c(
    1,
    3,
    5,
    10,
    NA
  )
)


test_that("filter_db filters by assembly level", {

  result <- filter_db(
    db,
    assembly_level = c(
      "Complete Genome",
      "Chromosome"
    ),
    report = FALSE
  )

  expect_equal(
    result$accession,
    c("A", "B", "E")
  )
})


test_that("filter_db filters by minimum contig N50", {

  result <- filter_db(
    db,
    min_contig_n50 = 50000,
    report = FALSE
  )

  expect_equal(
    result$accession,
    c("A", "B", "E")
  )
})


test_that("filter_db filters by minimum scaffold N50", {

  result <- filter_db(
    db,
    min_scaffold_n50 = 50000,
    report = FALSE
  )

  expect_equal(
    result$accession,
    c("A", "B", "C", "E")
  )
})


test_that("filter_db filters by minimum CheckM completeness", {

  result <- filter_db(
    db,
    min_checkm_completeness = 90,
    report = FALSE
  )

  expect_equal(
    result$accession,
    c("A", "B", "C", "E")
  )
})


test_that("filter_db filters by maximum CheckM contamination", {

  result <- filter_db(
    db,
    max_checkm_contamination = 5,
    report = FALSE
  )

  expect_equal(
    result$accession,
    c("A", "B", "C", "E")
  )
})


test_that("filter_db removes missing quality values when requested", {

  result <- filter_db(
    db,
    min_contig_n50 = 50000,
    keep_missing = FALSE,
    report = FALSE
  )

  expect_equal(
    result$accession,
    c("A", "B")
  )
})


test_that("filter_db combines multiple filters", {

  result <- filter_db(
    db,
    assembly_level = c(
      "Complete Genome",
      "Chromosome"
    ),
    min_contig_n50 = 50000,
    min_checkm_completeness = 95,
    max_checkm_contamination = 3,
    keep_missing = FALSE,
    report = FALSE
  )

  expect_equal(
    result$accession,
    c("A", "B")
  )
})


test_that("filter_db errors when required database columns are missing", {

  bad_db <- db |>
    dplyr::select(-accession)

  expect_error(
    filter_db(
      bad_db,
      report = FALSE
    ),
    "missing required column"
  )
})


test_that("filter_db errors when a requested quality column is missing", {

  bad_db <- db |>
    dplyr::select(-contig_n50)

  expect_error(
    filter_db(
      bad_db,
      min_contig_n50 = 50000,
      report = FALSE
    ),
    "Contig N50"
  )
})


test_that("filter_db validates keep_missing", {

  expect_error(
    filter_db(
      db,
      keep_missing = "yes",
      report = FALSE
    ),
    "keep_missing must be TRUE or FALSE"
  )
})


test_that("filter_db rejects invalid CheckM thresholds", {

  expect_error(
    filter_db(
      db,
      min_checkm_completeness = 150,
      report = FALSE
    )
  )

  expect_error(
    filter_db(
      db,
      max_checkm_contamination = -1,
      report = FALSE
    )
  )
})


test_that("filter_db leaves database unchanged when no filters are requested", {

  result <- filter_db(
    db,
    report = FALSE
  )

  expect_equal(
    result,
    db
  )
})
