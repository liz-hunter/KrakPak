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
  level = c(
    "S",
    "S",
    "S"
  ),
  db = c(
    "db1",
    "db2",
    "db1"
  ),
  perc_comp_exact = c(
    10,
    20,
    5
  ),
  incl_min_count = c(
    10,
    20,
    5
  )
)


# Presence ---------------------------------------------------------------

test_that("inspect_stats returns one row per taxon", {

  result <- inspect_stats(
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
})


test_that("inspect_stats reports replicate presence", {

  result <- inspect_stats(
    inspect_data,
    missing = "zero"
  )

  species_a <- result[
    result$taxid == 1001,
  ]

  species_b <- result[
    result$taxid == 1002,
  ]

  expect_equal(
    species_a$n_replicates,
    2
  )

  expect_equal(
    species_a$n_reported,
    2
  )

  expect_equal(
    species_a$report_rate,
    1
  )

  expect_equal(
    species_b$n_replicates,
    2
  )

  expect_equal(
    species_b$n_reported,
    1
  )

  expect_equal(
    species_b$report_rate,
    0.5
  )
})


# Missing taxa ------------------------------------------------------------

test_that("inspect_stats treats absent taxa as zero by default", {

  result <- inspect_stats(
    inspect_data,
    missing = "zero"
  )

  species_b <- result[
    result$taxid == 1002,
  ]

  # Species B is 5 in db1 and absent from db2,
  # so its values are treated as 5 and 0.
  expect_equal(
    species_b$mean,
    2.5
  )

  expect_equal(
    species_b$median,
    2.5
  )

  expect_equal(
    species_b$min,
    0
  )

  expect_equal(
    species_b$max,
    5
  )
})


test_that("inspect_stats can ignore absent taxa", {

  result <- inspect_stats(
    inspect_data,
    missing = "ignore"
  )

  species_b <- result[
    result$taxid == 1002,
  ]

  expect_equal(
    species_b$mean,
    5
  )

  expect_equal(
    species_b$n_reported,
    1
  )

  expect_equal(
    species_b$report_rate,
    0.5
  )
})


# Descriptive statistics -------------------------------------------------

test_that("inspect_stats calculates descriptive statistics correctly", {

  result <- inspect_stats(
    inspect_data
  )

  species_a <- result[
    result$taxid == 1001,
  ]

  expect_equal(
    species_a$mean,
    15
  )

  expect_equal(
    species_a$median,
    15
  )

  expect_equal(
    species_a$sd,
    stats::sd(c(10, 20))
  )

  expect_equal(
    species_a$min,
    10
  )

  expect_equal(
    species_a$max,
    20
  )

  expect_equal(
    species_a$range,
    10
  )

  expect_equal(
    species_a$iqr,
    stats::IQR(c(10, 20))
  )

  expect_equal(
    species_a$cv,
    100 * stats::sd(c(10, 20)) / 15
  )
})


# InfRV ------------------------------------------------------------------

test_that("inspect_stats calculates InfRV correctly", {

  result <- inspect_stats(
    inspect_data,
    pseudocount = 5
  )

  species_a <- result[
    result$taxid == 1001,
  ]

  # Counts are 10 and 20:
  # mean = 15
  # sample variance = 50
  # InfRV = (50 - 15) / (15 + 5) + 0.01
  expected <- (
    (50 - 15) /
      (15 + 5)
  ) + 0.01

  expect_equal(
    species_a$infrv_mean,
    15
  )

  expect_equal(
    species_a$infrv_variance,
    50
  )

  expect_equal(
    species_a$infrv,
    expected
  )
})


test_that("inspect_stats uses the requested pseudocount", {

  result <- inspect_stats(
    inspect_data,
    pseudocount = 10
  )

  species_a <- result[
    result$taxid == 1001,
  ]

  expected <- (
    (50 - 15) /
      (15 + 10)
  ) + 0.01

  expect_equal(
    species_a$infrv,
    expected
  )
})


# Validation --------------------------------------------------------------

test_that("inspect_stats requires expected columns", {

  bad_data <- inspect_data |>
    dplyr::select(-perc_comp_exact)

  expect_error(
    inspect_stats(
      bad_data
    ),
    "Missing required column"
  )
})


test_that("inspect_stats requires numeric value columns", {

  bad_data <- inspect_data

  bad_data$perc_comp_exact <- as.character(
    bad_data$perc_comp_exact
  )

  expect_error(
    inspect_stats(
      bad_data
    ),
    "value.*numeric column"
  )
})


test_that("inspect_stats rejects NA values", {

  bad_data <- inspect_data

  bad_data$perc_comp_exact[1] <- NA

  expect_error(
    inspect_stats(
      bad_data
    ),
    "contains NA values"
  )
})


test_that("inspect_stats rejects duplicate taxa within a replicate", {

  bad_data <- dplyr::bind_rows(
    inspect_data,
    inspect_data[1, ]
  )

  expect_error(
    inspect_stats(
      bad_data
    ),
    "more than once within a replicate"
  )
})


test_that("inspect_stats validates pseudocount", {

  expect_error(
    inspect_stats(
      inspect_data,
      pseudocount = -1
    ),
    "non-negative"
  )

  expect_error(
    inspect_stats(
      inspect_data,
      pseudocount = NA
    ),
    "non-negative"
  )
})
