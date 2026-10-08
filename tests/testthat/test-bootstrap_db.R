# Standard test data ------------------------------------------------------

db <- tibble::tibble(
  accession = c(
    "GCF_000001.1",
    "GCF_000002.1",
    "GCF_000003.1",
    "GCF_000004.1"
  ),
  name = c(
    "Genome A",
    "Genome B",
    "Genome C",
    "Genome D"
  ),
  taxid = c(
    101,
    102,
    102,
    103
  ),
  filename = c(
    "a.fna",
    "b.fna",
    "c.fna",
    "d.fna"
  )
)


# Basic behavior ----------------------------------------------------------

test_that("bootstrap_db returns expected run information", {

  outdir <- tempfile("bootstrap_")

  result <- bootstrap_db(
    db,
    n_boot = 3,
    label = "test",
    seed = 42,
    outdir = outdir,
    outputs = "composite",
    verbose = FALSE
  )

  expect_equal(result$label, "test")
  expect_equal(result$label_safe, "test")
  expect_equal(result$n_assemblies, 4)
  expect_equal(result$n_boot, 3)
  expect_equal(result$seed, 42)
  expect_equal(result$outputs, "composite")
})


test_that("bootstrap_db creates the requested output directory", {

  outdir <- tempfile("bootstrap_")

  bootstrap_db(
    db,
    n_boot = 2,
    label = "test",
    seed = 42,
    outdir = outdir,
    outputs = "composite",
    verbose = FALSE
  )

  expect_true(
    dir.exists(outdir)
  )
})


# Composite outputs -------------------------------------------------------

test_that("bootstrap_db writes composite accession and taxid matrices", {

  outdir <- tempfile("bootstrap_")

  result <- bootstrap_db(
    db,
    n_boot = 3,
    label = "test",
    seed = 42,
    outdir = outdir,
    outputs = "composite",
    verbose = FALSE
  )

  expect_true(
    file.exists(result$composite_accessions)
  )

  expect_true(
    file.exists(result$composite_taxids)
  )

  accessions <- readr::read_tsv(
    result$composite_accessions,
    show_col_types = FALSE
  )

  taxids <- readr::read_tsv(
    result$composite_taxids,
    show_col_types = FALSE
  )

  # Each bootstrap contains the same number of draws
  # as assemblies in the original database.
  expect_equal(
    nrow(accessions),
    nrow(db)
  )

  expect_equal(
    nrow(taxids),
    nrow(db)
  )

  # One column per bootstrap replicate.
  expect_equal(
    ncol(accessions),
    3
  )

  expect_equal(
    ncol(taxids),
    3
  )

  expect_equal(
    names(accessions),
    c(
      "test_db1",
      "test_db2",
      "test_db3"
    )
  )

  expect_equal(
    names(taxids),
    c(
      "test_db1",
      "test_db2",
      "test_db3"
    )
  )
})


test_that("bootstrap_db samples only assemblies from the input database", {

  outdir <- tempfile("bootstrap_")

  result <- bootstrap_db(
    db,
    n_boot = 5,
    label = "test",
    seed = 42,
    outdir = outdir,
    outputs = "composite",
    verbose = FALSE
  )

  accessions <- readr::read_tsv(
    result$composite_accessions,
    show_col_types = FALSE
  )

  expect_true(
    all(
      unlist(accessions) %in% db$accession
    )
  )
})


# Reproducibility ---------------------------------------------------------

test_that("bootstrap_db is reproducible when seed is specified", {

  outdir1 <- tempfile("bootstrap_")
  outdir2 <- tempfile("bootstrap_")

  result1 <- bootstrap_db(
    db,
    n_boot = 5,
    label = "test",
    seed = 123,
    outdir = outdir1,
    outputs = "composite",
    verbose = FALSE
  )

  result2 <- bootstrap_db(
    db,
    n_boot = 5,
    label = "test",
    seed = 123,
    outdir = outdir2,
    outputs = "composite",
    verbose = FALSE
  )

  accessions1 <- readr::read_tsv(
    result1$composite_accessions,
    show_col_types = FALSE
  )

  accessions2 <- readr::read_tsv(
    result2$composite_accessions,
    show_col_types = FALSE
  )

  expect_equal(
    accessions1,
    accessions2
  )
})


# Raw and unique files ----------------------------------------------------

test_that("bootstrap_db writes one raw file per bootstrap replicate", {

  outdir <- tempfile("bootstrap_")

  result <- bootstrap_db(
    db,
    n_boot = 4,
    label = "test",
    seed = 42,
    outdir = outdir,
    outputs = "files",
    verbose = FALSE
  )

  expect_length(
    result$bootstrap_files,
    4
  )

  expect_true(
    all(file.exists(result$bootstrap_files))
  )

  draws <- lapply(
    result$bootstrap_files,
    readLines
  )

  expect_true(
    all(
      lengths(draws) == nrow(db)
    )
  )
})


test_that("unique bootstrap files contain no duplicated accessions", {

  outdir <- tempfile("bootstrap_")

  result <- bootstrap_db(
    db,
    n_boot = 5,
    label = "test",
    seed = 42,
    outdir = outdir,
    outputs = "unique_files",
    verbose = FALSE
  )

  expect_length(
    result$unique_bootstrap_files,
    5
  )

  for (path in result$unique_bootstrap_files) {

    accessions <- readLines(path)

    expect_false(
      anyDuplicated(accessions) > 0
    )
  }
})


# Count tables ------------------------------------------------------------

test_that("accession count tables preserve total bootstrap draws", {

  outdir <- tempfile("bootstrap_")

  result <- bootstrap_db(
    db,
    n_boot = 3,
    label = "test",
    seed = 42,
    outdir = outdir,
    outputs = "accession_counts",
    verbose = FALSE
  )

  expect_length(
    result$accession_count_tables,
    3
  )

  for (path in result$accession_count_tables) {

    counts <- readr::read_tsv(
      path,
      show_col_types = FALSE
    )

    expect_equal(
      sum(counts$count),
      nrow(db)
    )
  }
})


test_that("taxid count tables preserve total bootstrap draws", {

  outdir <- tempfile("bootstrap_")

  result <- bootstrap_db(
    db,
    n_boot = 3,
    label = "test",
    seed = 42,
    outdir = outdir,
    outputs = "taxid_counts",
    verbose = FALSE
  )

  expect_length(
    result$taxid_count_tables,
    3
  )

  for (path in result$taxid_count_tables) {

    counts <- readr::read_tsv(
      path,
      show_col_types = FALSE
    )

    expect_equal(
      sum(counts$count),
      nrow(db)
    )
  }
})


# All-accessions output ---------------------------------------------------

test_that("all_accessions contains the union of sampled accessions", {

  outdir <- tempfile("bootstrap_")

  result <- bootstrap_db(
    db,
    n_boot = 5,
    label = "test",
    seed = 42,
    outdir = outdir,
    outputs = c(
      "composite",
      "all_accessions"
    ),
    verbose = FALSE
  )

  composite <- readr::read_tsv(
    result$composite_accessions,
    show_col_types = FALSE
  )

  expected <- sort(
    unique(unlist(composite))
  )

  expect_equal(
    result$all_accessions,
    expected
  )

  expect_true(
    file.exists(result$all_accessions_file)
  )

  expect_equal(
    readLines(result$all_accessions_file),
    expected
  )
})


# Label handling ----------------------------------------------------------

test_that("bootstrap_db sanitizes labels for output filenames", {

  outdir <- tempfile("bootstrap_")

  result <- bootstrap_db(
    db,
    n_boot = 2,
    label = "My test: database",
    seed = 42,
    outdir = outdir,
    outputs = "composite",
    verbose = FALSE
  )

  expect_equal(
    result$label_safe,
    "My_test__database"
  )

  expect_true(
    grepl(
      "My_test__database",
      basename(result$composite_accessions)
    )
  )
})


# Validation --------------------------------------------------------------

test_that("bootstrap_db validates n_boot", {

  expect_error(
    bootstrap_db(
      db,
      n_boot = 0,
      label = "test",
      outdir = tempfile(),
      verbose = FALSE
    ),
    "n_boot must be a positive integer"
  )

  expect_error(
    bootstrap_db(
      db,
      n_boot = 2.5,
      label = "test",
      outdir = tempfile(),
      verbose = FALSE
    ),
    "n_boot must be a positive integer"
  )
})


test_that("bootstrap_db validates output types", {

  expect_error(
    bootstrap_db(
      db,
      n_boot = 2,
      label = "test",
      outdir = tempfile(),
      outputs = "not_a_real_output",
      verbose = FALSE
    ),
    "Unknown output"
  )
})


test_that("bootstrap_db rejects duplicated input accessions", {

  bad_db <- dplyr::bind_rows(
    db,
    db[1, ]
  )

  expect_error(
    bootstrap_db(
      bad_db,
      n_boot = 2,
      label = "test",
      outdir = tempfile(),
      verbose = FALSE
    ),
    "duplicated accessions"
  )
})


test_that("bootstrap_db rejects missing required columns", {

  bad_db <- db |>
    dplyr::select(-filename)

  expect_error(
    bootstrap_db(
      bad_db,
      n_boot = 2,
      label = "test",
      outdir = tempfile(),
      verbose = FALSE
    ),
    "missing required column"
  )
})
