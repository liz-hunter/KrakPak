# Helper -----------------------------------------------------------------

write_test_inspect <- function(path, scale = 1) {

  lines <- c(
    paste(100, 100 * scale, 10 * scale, "R", 1, "root", sep = "\t"),
    paste(90,  90 * scale,  10 * scale, "D", 2, "  Bacteria", sep = "\t"),
    paste(80,  80 * scale,  20 * scale, "P", 1224, "    Proteobacteria", sep = "\t"),
    paste(60,  60 * scale,  10 * scale, "G", 1000, "      Genus alpha", sep = "\t"),
    paste(30,  30 * scale,  10 * scale, "S", 1001, "        Species alpha", sep = "\t"),
    paste(20,  20 * scale,  20 * scale, "S1", 1002, "          Strain alpha", sep = "\t"),
    paste(20,  20 * scale,  20 * scale, "S", 1003, "        Species beta", sep = "\t")
  )

  writeLines(lines, path)
}


make_inspect_directory <- function() {

  directory <- tempfile("inspect_")
  dir.create(directory)

  write_test_inspect(
    file.path(directory, "test_db1_inspect.txt"),
    scale = 1
  )

  write_test_inspect(
    file.path(directory, "test_db2_inspect.txt"),
    scale = 2
  )

  directory
}


# Basic reading -----------------------------------------------------------

test_that("read_inspect reads and combines matching reports", {

  directory <- make_inspect_directory()

  result <- read_inspect(
    directory,
    prefix = "test_",
    level = "species"
  )

  expect_equal(
    nrow(result),
    4
  )

  expect_equal(
    sort(unique(result$db)),
    c("test_db1", "test_db2")
  )

  expect_true(
    all(result$level == "S")
  )
})


test_that("read_inspect recalculates exact percentages", {

  directory <- make_inspect_directory()

  result <- read_inspect(
    directory,
    prefix = "test_",
    level = "species"
  )

  species_alpha <- result[
    result$taxid == 1001,
  ]

  expect_equal(
    species_alpha$perc_comp_exact,
    c(30, 30)
  )

  species_beta <- result[
    result$taxid == 1003,
  ]

  expect_equal(
    species_beta$perc_comp_exact,
    c(20, 20)
  )
})


test_that("read_inspect preserves report order and indentation", {

  directory <- make_inspect_directory()

  result <- read_inspect(
    directory,
    prefix = "test_db1",
    level = "all",
    fuzzy = TRUE
  )

  expect_equal(
    result$report_row,
    1:7
  )

  expect_equal(
    result$indent,
    c(0, 2, 4, 6, 8, 10, 8)
  )
})


# Rank filtering ----------------------------------------------------------

test_that("read_inspect excludes intermediate ranks when fuzzy is FALSE", {

  directory <- make_inspect_directory()

  result <- read_inspect(
    directory,
    prefix = "test_db1",
    level = "species",
    fuzzy = FALSE
  )

  expect_equal(
    result$taxid,
    c(1001, 1003)
  )

  expect_true(
    all(result$level == "S")
  )
})


test_that("read_inspect includes intermediate ranks when fuzzy is TRUE", {

  directory <- make_inspect_directory()

  result <- read_inspect(
    directory,
    prefix = "test_db1",
    level = "species",
    fuzzy = TRUE
  )

  expect_equal(
    result$taxid,
    c(1001, 1002, 1003)
  )

  expect_equal(
    result$level,
    c("S", "S1", "S")
  )
})


test_that("read_inspect returns all ranks when requested", {

  directory <- make_inspect_directory()

  result <- read_inspect(
    directory,
    prefix = "test_db1",
    level = "all",
    fuzzy = TRUE
  )

  expect_equal(
    nrow(result),
    7
  )

  expect_true(
    "S1" %in% result$level
  )
})


# Taxid inclusion and exclusion ------------------------------------------

test_that("read_inspect can restrict output to specified taxids", {

  directory <- make_inspect_directory()

  result <- read_inspect(
    directory,
    prefix = "test_",
    level = "all",
    fuzzy = TRUE,
    set_taxids = c(1001, 1003)
  )

  expect_equal(
    sort(unique(result$taxid)),
    c(1001, 1003)
  )
})


test_that("read_inspect accepts a data frame of taxids", {

  directory <- make_inspect_directory()

  taxids <- tibble::tibble(
    taxid = c(1001, 1003)
  )

  result <- read_inspect(
    directory,
    prefix = "test_",
    level = "all",
    fuzzy = TRUE,
    set_taxids = taxids
  )

  expect_equal(
    sort(unique(result$taxid)),
    c(1001, 1003)
  )
})


test_that("read_inspect can exclude specified taxids", {

  directory <- make_inspect_directory()

  result <- read_inspect(
    directory,
    prefix = "test_db1",
    level = "species",
    exclude_taxids = 1003
  )

  expect_equal(
    result$taxid,
    1001
  )
})


test_that("read_inspect rejects overlapping inclusion and exclusion taxids", {

  directory <- make_inspect_directory()

  expect_error(
    read_inspect(
      directory,
      prefix = "test_",
      level = "all",
      set_taxids = c(1001, 1003),
      exclude_taxids = 1003
    ),
    "present in both"
  )
})


# Validation --------------------------------------------------------------

test_that("read_inspect rejects an invalid directory", {

  expect_error(
    read_inspect(
      "directory_that_does_not_exist",
      prefix = "test_"
    ),
    "Directory does not exist"
  )
})


test_that("read_inspect errors when no files match the prefix", {

  directory <- make_inspect_directory()

  expect_error(
    read_inspect(
      directory,
      prefix = "not_a_real_prefix"
    ),
    "No files beginning with"
  )
})


test_that("read_inspect rejects invalid taxonomic levels", {

  directory <- make_inspect_directory()

  expect_error(
    read_inspect(
      directory,
      prefix = "test_",
      level = "superduperphylum"
    ),
    "level.*must be one of"
  )
})


test_that("read_inspect requires exactly one root row", {

  directory <- tempfile("inspect_")
  dir.create(directory)

  path <- file.path(
    directory,
    "bad_inspect.txt"
  )

  lines <- c(
    paste(100, 100, 10, "R", 1, "root", sep = "\t"),
    paste(100, 100, 10, "R", 1, "root again", sep = "\t")
  )

  writeLines(lines, path)

  expect_error(
    read_inspect(
      directory,
      prefix = "bad",
      level = "all"
    ),
    "exactly one root row"
  )
})
