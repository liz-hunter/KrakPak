# Standard test data ------------------------------------------------------

inspect_data <- tibble::tibble(
  taxid = c(
    1001,
    1001,
    1002
  ),
  name = c(
    "Species A",
    "Species A",
    "Species B"
  ),
  db = c(
    "db1",
    "db2",
    "db1"
  ),
  incl_min_count = c(
    10,
    20,
    5
  ),
  perc_comp_exact = c(
    10,
    20,
    5
  )
)


# Basic reshaping ---------------------------------------------------------

test_that("inspect_wide creates one row per taxon", {

  result <- inspect_wide(
    inspect_data
  )

  expect_equal(
    nrow(result),
    2
  )

  expect_equal(
    result$taxid,
    c(1001, 1002)
  )

  expect_equal(
    result$name,
    c("Species A", "Species B")
  )
})


test_that("inspect_wide creates one column per replicate", {

  result <- inspect_wide(
    inspect_data
  )

  expect_true(
    all(
      c(
        "taxid",
        "name",
        "db1",
        "db2"
      ) %in% names(result)
    )
  )
})


test_that("inspect_wide places replicate values correctly", {

  result <- inspect_wide(
    inspect_data
  )

  species_a <- result[
    result$taxid == 1001,
  ]

  expect_equal(
    species_a$db1,
    10
  )

  expect_equal(
    species_a$db2,
    20
  )
})


# Missing combinations ---------------------------------------------------

test_that("inspect_wide fills absent taxa with zero by default", {

  result <- inspect_wide(
    inspect_data
  )

  species_b <- result[
    result$taxid == 1002,
  ]

  expect_equal(
    species_b$db1,
    5
  )

  expect_equal(
    species_b$db2,
    0
  )
})


test_that("inspect_wide can retain absent combinations as NA", {

  result <- inspect_wide(
    inspect_data,
    fill = NA
  )

  species_b <- result[
    result$taxid == 1002,
  ]

  expect_true(
    is.na(species_b$db2)
  )
})


# Alternative value columns ----------------------------------------------

test_that("inspect_wide can reshape a different value column", {

  result <- inspect_wide(
    inspect_data,
    value = "perc_comp_exact"
  )

  species_a <- result[
    result$taxid == 1001,
  ]

  expect_equal(
    species_a$db1,
    10
  )

  expect_equal(
    species_a$db2,
    20
  )
})


# Validation --------------------------------------------------------------

test_that("inspect_wide requires expected columns", {

  bad_data <- inspect_data |>
    dplyr::select(-incl_min_count)

  expect_error(
    inspect_wide(
      bad_data
    ),
    "Missing required column"
  )
})


test_that("inspect_wide requires a numeric value column", {

  bad_data <- inspect_data

  bad_data$incl_min_count <- as.character(
    bad_data$incl_min_count
  )

  expect_error(
    inspect_wide(
      bad_data
    ),
    "numeric column"
  )
})


test_that("inspect_wide rejects duplicate taxa within a replicate", {

  bad_data <- dplyr::bind_rows(
    inspect_data,
    inspect_data[1, ]
  )

  expect_error(
    inspect_wide(
      bad_data
    ),
    "more than once within a replicate"
  )
})
