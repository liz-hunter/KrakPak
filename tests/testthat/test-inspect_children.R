# Helper -----------------------------------------------------------------

make_children_inspect <- function() {

  directory <- tempfile("inspect_children_")
  dir.create(directory)

  lines <- c(
    paste(100, 100, 10, "R", 1, "root", sep = "\t"),
    paste(90,  90, 10, "D", 2, "  Bacteria", sep = "\t"),
    paste(80,  80, 20, "P", 1224, "    Proteobacteria", sep = "\t"),
    paste(60,  60, 10, "G", 1000, "      Genus alpha", sep = "\t"),
    paste(30,  30, 10, "S", 1001, "        Species alpha", sep = "\t"),
    paste(20,  20, 20, "S1", 1002, "          Strain alpha", sep = "\t"),
    paste(20,  20, 20, "S", 1003, "        Species beta", sep = "\t")
  )

  writeLines(
    lines,
    file.path(
      directory,
      "test_db1_inspect.txt"
    )
  )

  read_inspect(
    directory,
    prefix = "test_",
    level = "all",
    fuzzy = TRUE
  )
}


# Descendant extraction ---------------------------------------------------

test_that("inspect_children returns all descendants", {

  inspect <- make_children_inspect()

  result <- inspect_children(
    inspect,
    taxid = 1000
  )

  expect_equal(
    result$taxid,
    c(
      1000,
      1001,
      1002,
      1003
    )
  )
})


test_that("inspect_children can exclude the parent", {

  inspect <- make_children_inspect()

  result <- inspect_children(
    inspect,
    taxid = 1000,
    include_parent = FALSE
  )

  expect_equal(
    result$taxid,
    c(
      1001,
      1002,
      1003
    )
  )
})


test_that("inspect_children returns only immediate children when requested", {

  inspect <- make_children_inspect()

  result <- inspect_children(
    inspect,
    taxid = 1000,
    immediate = TRUE
  )

  expect_equal(
    result$taxid,
    c(
      1000,
      1001,
      1003
    )
  )
})


test_that("inspect_children identifies intermediate-rank children correctly", {

  inspect <- make_children_inspect()

  result <- inspect_children(
    inspect,
    taxid = 1001,
    immediate = TRUE
  )

  expect_equal(
    result$taxid,
    c(
      1001,
      1002
    )
  )
})


# Tree validation ---------------------------------------------------------

test_that("inspect_children validates minimizer counts", {

  inspect <- make_children_inspect()

  expect_no_error(
    inspect_children(
      inspect,
      taxid = 1000,
      validate = TRUE
    )
  )
})


test_that("inspect_children detects an inconsistent hierarchy", {

  inspect <- make_children_inspect()

  inspect$excl_min_count[
    inspect$taxid == 1003
  ] <- 999

  expect_error(
    inspect_children(
      inspect,
      taxid = 1000,
      validate = TRUE
    ),
    "Tree validation failed"
  )
})


test_that("inspect_children can skip tree validation", {

  inspect <- make_children_inspect()

  inspect$excl_min_count[
    inspect$taxid == 1003
  ] <- 999

  expect_no_error(
    inspect_children(
      inspect,
      taxid = 1000,
      validate = FALSE
    )
  )
})


# Validation --------------------------------------------------------------

test_that("inspect_children errors when taxid is absent", {

  inspect <- make_children_inspect()

  expect_error(
    inspect_children(
      inspect,
      taxid = 999999
    ),
    "was not found"
  )
})


test_that("inspect_children requires a single taxid", {

  inspect <- make_children_inspect()

  expect_error(
    inspect_children(
      inspect,
      taxid = c(1000, 1001)
    ),
    "exactly one"
  )
})


test_that("inspect_children requires hierarchy columns", {

  inspect <- make_children_inspect() |>
    dplyr::select(-indent)

  expect_error(
    inspect_children(
      inspect,
      taxid = 1000
    ),
    "missing required columns"
  )
})


test_that("inspect_children validates logical arguments", {

  inspect <- make_children_inspect()

  expect_error(
    inspect_children(
      inspect,
      taxid = 1000,
      include_parent = "yes"
    ),
    "include_parent"
  )

  expect_error(
    inspect_children(
      inspect,
      taxid = 1000,
      immediate = NA
    ),
    "immediate"
  )

  expect_error(
    inspect_children(
      inspect,
      taxid = 1000,
      validate = 1
    ),
    "validate"
  )
})
