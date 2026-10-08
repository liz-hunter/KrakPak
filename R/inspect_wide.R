#' Reshape Kraken inspect replicates to wide format
#'
#' Converts replicate Kraken inspect reports from long format to a wide
#' taxon-by-replicate table. Each row represents one taxon and each replicate
#' database becomes a separate column.
#'
#' This is useful for directly examining replicate-level values or exporting
#' matrices for downstream analyses.
#'
#' @param data A tibble returned by `read_inspect()`.
#' @param name The dataset name
#' @param value Numeric column to place in the replicate columns. Defaults to
#'   `"incl_min_count"`. Common alternatives include `"perc_comp_exact"` and
#'   `"excl_min_count"`.
#' @param replicate Column identifying replicate databases. Defaults to `"db"`.
#' @param taxid Column identifying taxa. Defaults to `"taxid"`.
#' @param fill Value used when a taxon is absent from a replicate. Defaults
#'   to `0`. Use `NA` to retain absent combinations as missing.
#'
#' @return A tibble with one row per taxon and one column per replicate.
#'
#' @export
inspect_wide <- function(
    data,
    value = "incl_min_count",
    replicate = "db",
    taxid = "taxid",
    name = "name",
    fill = 0
) {

  required_cols <- c(taxid, replicate, value)
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

  # Check for more than one value for a taxon within a replicate
  duplicates <- data |>
    dplyr::count(
      dplyr::across(
        dplyr::all_of(c(taxid, replicate))
      )
    ) |>
    dplyr::filter(n > 1)

  if (nrow(duplicates) > 0) {
    stop(
      "Some taxa occur more than once within a replicate.",
      call. = FALSE
    )
  }

  data |>
    dplyr::select(
      dplyr::all_of(c(taxid, name, replicate, value))
    ) |>
    tidyr::pivot_wider(
      names_from = dplyr::all_of(replicate),
      values_from = dplyr::all_of(value),
      values_fill = fill
    ) |>
    dplyr::arrange(.data[[taxid]])
}
