# Helper -----------------------------------------------------------------

make_test_boot <- function() {

  outdir <- tempfile("bootstrap_stats_")
  dir.create(outdir)

  taxids <- tibble::tibble(
    test_db1 = c(1, 1, 2, 2),
    test_db2 = c(1, 1, 1, 2),
    test_db3 = c(1, 1, 1, 1)
  )

  taxid_path <- file.path(
    outdir,
    "composite_test_taxids.tsv"
  )

  readr::write_tsv(
    taxids,
    taxid_path
  )

  list(
    label = "test",
    label_safe = "test",
    outdir = outdir,
    composite_taxids = taxid_path
  )
}


# Basic output ------------------------------------------------------------

test_that("bootstrap_stats returns one row per bootstrap replicate", {

  boot <- make_test_boot()

  result <- bootstrap_stats(
    boot,
    write = FALSE,
    verbose = FALSE
  )

  expect_equal(
    nrow(result$summary),
    3
  )

  expect_equal(
    result$summary$replicate,
    c(1, 2, 3)
  )

  expect_equal(
    result$summary$sample_id,
    c(
      "test_db1",
      "test_db2",
      "test_db3"
    )
  )
})


test_that("bootstrap_stats reports the expected columns", {

  boot <- make_test_boot()

  result <- bootstrap_stats(
    boot,
    write = FALSE,
    verbose = FALSE
  )

  expect_true(
    all(
      c(
        "sample_id",
        "replicate",
        "n_draws",
        "n_unique_taxa",
        "shannon",
        "simpson_dominance",
        "gini_simpson",
        "inverse_simpson",
        "pielou_evenness"
      ) %in% names(result$summary)
    )
  )
})


# Number of draws and taxa ------------------------------------------------

test_that("bootstrap_stats counts draws correctly", {

  boot <- make_test_boot()

  result <- bootstrap_stats(
    boot,
    write = FALSE,
    verbose = FALSE
  )

  expect_equal(
    result$summary$n_draws,
    c(4, 4, 4)
  )
})


test_that("bootstrap_stats counts unique taxa correctly", {

  boot <- make_test_boot()

  result <- bootstrap_stats(
    boot,
    write = FALSE,
    verbose = FALSE
  )

  expect_equal(
    result$summary$n_unique_taxa,
    c(2, 2, 1)
  )
})


# Diversity metrics -------------------------------------------------------

test_that("bootstrap_stats calculates Shannon diversity correctly", {

  boot <- make_test_boot()

  result <- bootstrap_stats(
    boot,
    write = FALSE,
    verbose = FALSE
  )

  # test_db1 has proportions 0.5 / 0.5.
  expected_db1 <- log(2)

  # test_db2 has proportions 0.75 / 0.25.
  expected_db2 <- -(
    0.75 * log(0.75) +
      0.25 * log(0.25)
  )

  # test_db3 contains only one taxon.
  expected_db3 <- 0

  expect_equal(
    result$summary$shannon,
    c(
      expected_db1,
      expected_db2,
      expected_db3
    ),
    tolerance = 1e-10
  )
})


test_that("bootstrap_stats calculates Simpson dominance correctly", {

  boot <- make_test_boot()

  result <- bootstrap_stats(
    boot,
    write = FALSE,
    verbose = FALSE
  )

  expect_equal(
    result$summary$simpson_dominance,
    c(
      0.5,
      0.625,
      1
    ),
    tolerance = 1e-10
  )
})


test_that("bootstrap_stats calculates Gini-Simpson diversity correctly", {

  boot <- make_test_boot()

  result <- bootstrap_stats(
    boot,
    write = FALSE,
    verbose = FALSE
  )

  expect_equal(
    result$summary$gini_simpson,
    c(
      0.5,
      0.375,
      0
    ),
    tolerance = 1e-10
  )
})


test_that("bootstrap_stats calculates inverse Simpson diversity correctly", {

  boot <- make_test_boot()

  result <- bootstrap_stats(
    boot,
    write = FALSE,
    verbose = FALSE
  )

  expect_equal(
    result$summary$inverse_simpson,
    c(
      2,
      1.6,
      1
    ),
    tolerance = 1e-10
  )
})


test_that("bootstrap_stats calculates Pielou evenness correctly", {

  boot <- make_test_boot()

  result <- bootstrap_stats(
    boot,
    write = FALSE,
    verbose = FALSE
  )

  expected_db2 <- -(
    0.75 * log(0.75) +
      0.25 * log(0.25)
  ) / log(2)

  expect_equal(
    result$summary$pielou_evenness[1],
    1,
    tolerance = 1e-10
  )

  expect_equal(
    result$summary$pielou_evenness[2],
    expected_db2,
    tolerance = 1e-10
  )

  # Evenness is undefined when there is only one taxon.
  expect_true(
    is.na(
      result$summary$pielou_evenness[3]
    )
  )
})


# File output -------------------------------------------------------------

test_that("bootstrap_stats does not write a file when write is FALSE", {

  boot <- make_test_boot()

  result <- bootstrap_stats(
    boot,
    write = FALSE,
    verbose = FALSE
  )

  expect_null(
    result$path
  )
})


test_that("bootstrap_stats writes the summary when requested", {

  boot <- make_test_boot()

  metrics_dir <- tempfile("metrics_")

  result <- bootstrap_stats(
    boot,
    outdir = metrics_dir,
    write = TRUE,
    verbose = FALSE
  )

  expect_true(
    file.exists(result$path)
  )

  written <- readr::read_tsv(
    result$path,
    show_col_types = FALSE
  )

  expect_equal(
    written,
    result$summary
  )
})


# Integration with bootstrap_db ------------------------------------------

test_that("bootstrap_stats works with bootstrap_db output", {

  db <- tibble::tibble(
    accession = c(
      "A",
      "B",
      "C",
      "D"
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

  boot <- bootstrap_db(
    db,
    n_boot = 5,
    label = "integration",
    seed = 42,
    outdir = tempfile("bootstrap_"),
    outputs = "composite",
    verbose = FALSE
  )

  result <- bootstrap_stats(
    boot,
    write = FALSE,
    verbose = FALSE
  )

  expect_equal(
    nrow(result$summary),
    5
  )

  expect_true(
    all(
      result$summary$n_draws == nrow(db)
    )
  )
})


# Validation --------------------------------------------------------------

test_that("bootstrap_stats requires a bootstrap result", {

  expect_error(
    bootstrap_stats(
      "not a bootstrap result",
      write = FALSE,
      verbose = FALSE
    ),
    "boot must be the result returned by bootstrap_db"
  )
})


test_that("bootstrap_stats requires composite taxid output", {

  bad_boot <- list(
    label = "test",
    label_safe = "test",
    outdir = tempfile()
  )

  expect_error(
    bootstrap_stats(
      bad_boot,
      write = FALSE,
      verbose = FALSE
    ),
    "No composite taxid file found"
  )
})


test_that("bootstrap_stats detects a missing composite taxid file", {

  bad_boot <- list(
    label = "test",
    label_safe = "test",
    outdir = tempfile(),
    composite_taxids = "this_file_does_not_exist.tsv"
  )

  expect_error(
    bootstrap_stats(
      bad_boot,
      write = FALSE,
      verbose = FALSE
    ),
    "Composite taxid file does not exist"
  )
})
