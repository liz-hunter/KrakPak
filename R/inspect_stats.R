#' Calculate variability statistics across Kraken inspect replicates
#'
#' Summarizes taxon representation across replicate Kraken database inspect
#' reports and calculates measures of replicate variability.
#'
#' Inferential relative variance (InfRV) is calculated following Zhu et al.
#' (2019), using a count-valued column:
#'
#' \deqn{
#'   InfRV = \frac{\max(s^2 - \mu, 0)}{\mu + pseudocount} + 0.01
#' }
#'
#' where \eqn{s^2} is the sample variance and \eqn{\mu} is the mean across
#' replicate databases. The original publication used a pseudocount of 5.
#'
#' Taxa absent from individual replicate reports are treated as zero by
#' default.
#'
#' @param data A tibble returned by `read_inspect()`.
#' @param value Numeric column used for the standard summary statistics.
#'   Defaults to `"perc_comp_exact"`.
#' @param infrv_value Count-valued column used to calculate inferential
#'   relative variance. Defaults to `"incl_min_count"`.
#' @param pseudocount Non-negative numeric pseudocount added to the InfRV
#'   denominator. Defaults to 5, as used by Zhu et al. (2019).
#' @param replicate Column identifying replicate databases. Defaults to `"db"`.
#' @param taxon_cols Columns identifying taxa. Defaults to
#'   `c("taxid", "name", "level")`.
#' @param missing How taxa absent from individual replicates should be handled.
#'   `"zero"` treats absence as zero; `"ignore"` calculates statistics using
#'   only replicates in which the taxon appears.
#'
#' @return A tibble with one row per taxon containing presence, variability,
#'   and inferential relative variance statistics.
#'
#' @export
inspect_stats <- function(
    data,
    value = "perc_comp_exact",
    infrv_value = "incl_min_count",
    pseudocount = 5,
    replicate = "db",
    taxon_cols = c("taxid", "name", "level"),
    missing = c("zero", "ignore")
) {

  missing <- match.arg(missing)

  required_cols <- unique(
    c(value, infrv_value, replicate, taxon_cols)
  )

  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(data[[value]])) {
    stop(
      "`value` must identify a numeric column.",
      call. = FALSE
    )
  }

  if (!is.numeric(data[[infrv_value]])) {
    stop(
      "`infrv_value` must identify a numeric column.",
      call. = FALSE
    )
  }

  if (
    length(pseudocount) != 1 ||
    !is.numeric(pseudocount) ||
    is.na(pseudocount) ||
    pseudocount < 0
  ) {
    stop(
      "`pseudocount` must be a single non-negative number.",
      call. = FALSE
    )
  }

  x <- data |>
    dplyr::select(dplyr::all_of(required_cols)) |>
    dplyr::rename(
      .value = dplyr::all_of(value),
      .infrv_value = dplyr::all_of(infrv_value),
      .replicate = dplyr::all_of(replicate)
    ) |>
    dplyr::mutate(.reported = TRUE)

  if (anyNA(x$.value)) {
    stop(
      "`", value, "` contains NA values.",
      call. = FALSE
    )
  }

  if (anyNA(x$.infrv_value)) {
    stop(
      "`", infrv_value, "` contains NA values.",
      call. = FALSE
    )
  }

  n_replicates <- dplyr::n_distinct(x$.replicate)

  # Check for duplicate taxon observations within replicates
  duplicates <- x |>
    dplyr::count(
      dplyr::across(
        dplyr::all_of(c(taxon_cols, ".replicate"))
      )
    ) |>
    dplyr::filter(n > 1)

  if (nrow(duplicates) > 0) {
    stop(
      "Some taxa occur more than once within a replicate.",
      call. = FALSE
    )
  }

  # Add absent taxon/replicate combinations as zero
  if (missing == "zero") {

    taxa <- x |>
      dplyr::distinct(
        dplyr::across(dplyr::all_of(taxon_cols))
      )

    replicates <- x |>
      dplyr::distinct(.replicate)

    x <- tidyr::crossing(taxa, replicates) |>
      dplyr::left_join(
        x,
        by = c(taxon_cols, ".replicate")
      ) |>
      dplyr::mutate(
        .value = tidyr::replace_na(.value, 0),
        .infrv_value = tidyr::replace_na(.infrv_value, 0),
        .reported = tidyr::replace_na(.reported, FALSE)
      )
  }

  x |>
    dplyr::group_by(
      dplyr::across(dplyr::all_of(taxon_cols))
    ) |>
    dplyr::summarise(

      # Replicate presence
      n_replicates = n_replicates,
      n_reported = sum(.reported),
      report_rate = n_reported / n_replicates,

      # Descriptive statistics for `value`
      mean = mean(.value),
      median = stats::median(.value),
      sd = stats::sd(.value),

      cv = dplyr::if_else(
        mean > 0,
        100 * sd / mean,
        NA_real_
      ),

      min = min(.value),
      max = max(.value),
      range = max - min,
      iqr = stats::IQR(.value),

      # Components of InfRV
      infrv_mean = mean(.infrv_value),
      infrv_variance = stats::var(.infrv_value),

      # Zhu et al. (2019) InfRV
      infrv = (
        pmax(
          infrv_variance - infrv_mean,
          0
        ) /
          (infrv_mean + pseudocount)
      ) + 0.01,

      .groups = "drop"
    )
}
