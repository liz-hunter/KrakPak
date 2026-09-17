#' Calculate variability statistics across Kraken inspect replicates
#'
#' Calculates summary and variability statistics for taxa across replicate
#' Kraken database inspect reports. Useful for identifying taxa that are
#' consistently represented across databases versus taxa whose representation
#' is highly variable.
#'
#' @param data A tibble returned by `read_inspect()`.
#' @param value Numeric column to summarize. Defaults to `"perc_comp_exact"`.
#' @param replicate Column identifying replicate databases. Defaults to `"db"`.
#' @param taxon_cols Columns identifying taxa. Defaults to
#'   `c("taxid", "name", "level")`.
#' @param missing How taxa absent from individual replicates should be handled.
#'   `"zero"` treats absence as zero; `"ignore"` calculates statistics using
#'   only replicates in which the taxon appears.
#'
#' @return A tibble containing one row per taxon with measures of central
#'   tendency, dispersion, and replicate consistency.
#'
#' @export
inspect_stats <- function(
    data,
    value = "perc_comp_exact",
    replicate = "db",
    taxon_cols = c("taxid", "name", "level"),
    missing = c("zero", "ignore")
) {

  missing <- match.arg(missing)

  required_cols <- c(value, replicate, taxon_cols)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(data[[value]])) {
    stop("`value` must identify a numeric column.", call. = FALSE)
  }

  x <- data %>%
    dplyr::select(dplyr::all_of(required_cols)) %>%
    dplyr::rename(
      .value = dplyr::all_of(value),
      .replicate = dplyr::all_of(replicate)
    ) %>%
    dplyr::mutate(.reported = TRUE)

  if (anyNA(x$.value)) {
    stop(
      "`", value, "` contains NA values.",
      call. = FALSE
    )
  }

  n_replicates <- dplyr::n_distinct(x$.replicate)

  # Check that each taxon occurs at most once per replicate.
  duplicates <- x %>%
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(c(taxon_cols, ".replicate"))
      )
    ) %>%
    dplyr::summarise(
      .n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(.n > 1)

  if (nrow(duplicates) > 0) {
    stop(
      "Some taxa occur more than once within a replicate.",
      call. = FALSE
    )
  }

  # Add zeroes for taxa missing from individual replicates.
  if (missing == "zero") {

    taxa <- x %>%
      dplyr::distinct(
        dplyr::across(dplyr::all_of(taxon_cols))
      )

    replicates <- x %>%
      dplyr::distinct(.replicate)

    x <- tidyr::crossing(taxa, replicates) %>%
      dplyr::left_join(
        x,
        by = c(taxon_cols, ".replicate")
      ) %>%
      dplyr::mutate(
        .value = tidyr::replace_na(.value, 0),
        .reported = tidyr::replace_na(.reported, FALSE)
      )
  }

  x %>%
    dplyr::group_by(
      dplyr::across(dplyr::all_of(taxon_cols))
    ) %>%
    dplyr::summarise(

      # Replicate consistency
      n_replicates = n_replicates,
      n_reported = sum(.reported),
      report_rate = n_reported / n_replicates,

      n_nonzero = sum(.value > 0),
      nonzero_rate = n_nonzero / n_replicates,

      # Central tendency
      mean = mean(.value),
      median = stats::median(.value),

      # Classical dispersion
      sd = stats::sd(.value),
      variance = stats::var(.value),

      min = min(.value),
      max = max(.value),
      range = max - min,

      # Robust dispersion
      q1 = as.numeric(
        stats::quantile(.value, 0.25, names = FALSE)
      ),

      q3 = as.numeric(
        stats::quantile(.value, 0.75, names = FALSE)
      ),

      iqr = stats::IQR(.value),

      # Unscaled median absolute deviation
      mad = stats::mad(
        .value,
        constant = 1
      ),

      .groups = "drop"
    ) %>%
    dplyr::mutate(

      # Coefficient of variation
      cv = dplyr::if_else(
        mean != 0 & !is.na(sd),
        sd / abs(mean),
        NA_real_
      ),

      cv_percent = cv * 100,

      # Robust CV based on MAD.
      # 1.4826 scales MAD to approximate SD under normality.
      robust_cv = dplyr::if_else(
        median != 0,
        (1.4826 * mad) / abs(median),
        NA_real_
      ),

      robust_cv_percent = robust_cv * 100,

      # Total observed spread relative to the mean
      relative_range = dplyr::if_else(
        mean != 0,
        range / abs(mean),
        NA_real_
      ),

      relative_range_percent = relative_range * 100,

      # Robust relative dispersion
      quartile_dispersion = dplyr::if_else(
        (q1 + q3) != 0,
        (q3 - q1) / (q3 + q1),
        NA_real_
      ),

      # 1 = perfectly consistent min/max
      # 0 = absent/zero in at least one replicate
      min_max_ratio = dplyr::if_else(
        max > 0,
        min / max,
        NA_real_
      )
    )
}
