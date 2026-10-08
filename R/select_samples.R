#' Classify and select taxa for in silico evaluation
#'
#' Convenience wrapper around `classify_inspect_taxa()` and
#' `select_inspect_taxa()`.
#'
#' By default, low and high variability taxa are selected by ranking,
#' favoring better represented taxa, while intermediate variability taxa
#' are sampled randomly using a reproducible seed.
#'
#' @param stats A data frame produced by `inspect_stats()`.
#' @param tax_level Taxonomic level to retain. Defaults to `"S"`.
#' @param n_per_group Number of taxa to select from each group. May be a
#'   single integer or a named vector with entries `low`, `intermediate`,
#'   and `high`.
#' @param seed Integer random seed. Defaults to `20260917`.
#' @param low_report_min Minimum report rate for the low variability class.
#' @param low_cv_max Maximum coefficient of variation for the low
#'   variability class.
#' @param mid_report_min Minimum report rate for the intermediate class.
#' @param mid_report_max Maximum report rate for the intermediate class.
#' @param mid_cv_min Minimum coefficient of variation for the intermediate
#'   class.
#' @param mid_cv_max Maximum coefficient of variation for the intermediate
#'   class.
#' @param high_report_max Maximum report rate for the high variability class.
#' @param high_cv_min Minimum coefficient of variation for the high
#'   variability class.
#' @param min_mean Optional minimum mean representation required for all
#'   classes. Defaults to 0.
#' @param low_method Selection method for the low class: `"ranked"` or
#'   `"random"`.
#' @param intermediate_method Selection method for the intermediate class:
#'   `"ranked"` or `"random"`.
#' @param high_method Selection method for the high class: `"ranked"` or
#'   `"random"`.
#' @param strict Logical. If `TRUE`, error when fewer taxa are available
#'   than requested. If `FALSE`, return all available taxa with a warning.
#'
#' @return A list containing:
#'
#' * `classified`: all taxa at the requested taxonomic level with assigned
#'   variability classes.
#' * `candidates`: taxa satisfying low, intermediate, or high criteria.
#' * `selected`: taxa chosen for in silico evaluation.
#' * `summary`: number of candidate and selected taxa in each class.
#'
#' @export
classify_inspect_taxa <- function(
    stats,
    tax_level = "S",

    low_report_min = 0.9,
    low_cv_max = 40,

    mid_report_min = 0.6,
    mid_report_max = 0.8,
    mid_cv_min = 50,
    mid_cv_max = 100,

    high_report_max = 0.3,
    high_cv_min = 160,

    min_mean = 0
) {

  if (!is.data.frame(stats)) {
    stop(
      "stats must be a data frame produced by inspect_stats().",
      call. = FALSE
    )
  }

  required_cols <- c(
    "taxid",
    "name",
    "level",
    "n_replicates",
    "n_reported",
    "report_rate",
    "mean",
    "median",
    "sd",
    "cv",
    "min",
    "max",
    "range",
    "iqr"
  )

  missing_cols <- setdiff(
    required_cols,
    names(stats)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (
    length(tax_level) != 1L ||
    !is.character(tax_level) ||
    is.na(tax_level) ||
    !nzchar(tax_level)
  ) {
    stop(
      "tax_level must be a single non-empty character value.",
      call. = FALSE
    )
  }

  validate_report_rate <- function(x, arg) {

    if (
      length(x) != 1L ||
      !is.numeric(x) ||
      is.na(x) ||
      !is.finite(x) ||
      x < 0 ||
      x > 1
    ) {
      stop(
        arg,
        " must be a single numeric value between 0 and 1.",
        call. = FALSE
      )
    }

    x
  }

  validate_nonnegative <- function(x, arg) {

    if (
      length(x) != 1L ||
      !is.numeric(x) ||
      is.na(x) ||
      !is.finite(x) ||
      x < 0
    ) {
      stop(
        arg,
        " must be a single non-negative numeric value.",
        call. = FALSE
      )
    }

    x
  }

  low_report_min <- validate_report_rate(
    low_report_min,
    "low_report_min"
  )

  mid_report_min <- validate_report_rate(
    mid_report_min,
    "mid_report_min"
  )

  mid_report_max <- validate_report_rate(
    mid_report_max,
    "mid_report_max"
  )

  high_report_max <- validate_report_rate(
    high_report_max,
    "high_report_max"
  )

  low_cv_max <- validate_nonnegative(
    low_cv_max,
    "low_cv_max"
  )

  mid_cv_min <- validate_nonnegative(
    mid_cv_min,
    "mid_cv_min"
  )

  mid_cv_max <- validate_nonnegative(
    mid_cv_max,
    "mid_cv_max"
  )

  high_cv_min <- validate_nonnegative(
    high_cv_min,
    "high_cv_min"
  )

  min_mean <- validate_nonnegative(
    min_mean,
    "min_mean"
  )

  if (mid_report_min > mid_report_max) {
    stop(
      "mid_report_min must be <= mid_report_max.",
      call. = FALSE
    )
  }

  if (mid_cv_min > mid_cv_max) {
    stop(
      "mid_cv_min must be <= mid_cv_max.",
      call. = FALSE
    )
  }

  out <- tibble::as_tibble(stats)

  out <- out |>
    dplyr::filter(
      .data$level == tax_level
    ) |>
    dplyr::mutate(
      variability_class = dplyr::case_when(

        is.finite(.data$mean) &
          is.finite(.data$cv) &
          .data$mean >= min_mean &
          .data$report_rate >= low_report_min &
          .data$cv <= low_cv_max ~ "low",

        is.finite(.data$mean) &
          is.finite(.data$cv) &
          .data$mean >= min_mean &
          .data$report_rate >= mid_report_min &
          .data$report_rate <= mid_report_max &
          .data$cv >= mid_cv_min &
          .data$cv <= mid_cv_max ~ "intermediate",

        is.finite(.data$mean) &
          is.finite(.data$cv) &
          .data$mean >= min_mean &
          .data$report_rate <= high_report_max &
          .data$cv >= high_cv_min ~ "high",

        TRUE ~ "other"
      ),
      variability_class = factor(
        .data$variability_class,
        levels = c(
          "low",
          "intermediate",
          "high",
          "other"
        )
      )
    )

  out
}


#' Select taxa from variability classes
#'
#' Selects taxa from low, intermediate, and high variability classes
#' produced by classify_inspect_taxa().
#'
#' Each variability class may be selected either by ranking or by random
#' sampling. Ranked selection favors better represented taxa by sorting
#' primarily on mean representation.
#'
#' Random sampling is reproducible when a seed is supplied. Each
#' variability class receives an independent seed derived from the supplied
#' seed so that changing the selection method for one class does not alter
#' random selections in another class.
#'
#' @param classified A data frame produced by [classify_inspect_taxa()].
#' @param n_per_group Number of taxa to select from each group. May be a
#'   single integer or a named vector with entries `low`, `intermediate`,
#'   and `high`.
#' @param seed Integer random seed. Defaults to `20260917`.
#' @param low_method Selection method for the low class: `"ranked"` or
#'   `"random"`.
#' @param intermediate_method Selection method for the intermediate class:
#'   `"ranked"` or `"random"`.
#' @param high_method Selection method for the high class: `"ranked"` or
#'   `"random"`.
#' @param strict Logical. If `TRUE`, error when fewer taxa are available
#'   than requested. If `FALSE`, return all available taxa with a warning.
#'
#' @return A tibble containing the selected taxa plus `selection_method`
#'   and `selection_rank`.
#'
#' @export
select_inspect_taxa <- function(
    classified,
    n_per_group = 10,
    seed = 20260917L,
    low_method = "ranked",
    intermediate_method = "random",
    high_method = "ranked",
    strict = FALSE
) {

  if (!is.data.frame(classified)) {
    stop(
      "classified must be a data frame produced by classify_inspect_taxa().",
      call. = FALSE
    )
  }

  required_cols <- c(
    "taxid",
    "name",
    "mean",
    "cv",
    "report_rate",
    "variability_class"
  )

  missing_cols <- setdiff(
    required_cols,
    names(classified)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  groups <- c(
    "low",
    "intermediate",
    "high"
  )

  if (length(n_per_group) == 1L) {

    if (
      !is.numeric(n_per_group) ||
      is.na(n_per_group) ||
      !is.finite(n_per_group) ||
      n_per_group < 0 ||
      n_per_group != as.integer(n_per_group)
    ) {
      stop(
        "n_per_group must contain non-negative integers.",
        call. = FALSE
      )
    }

    n_per_group <- stats::setNames(
      rep(
        as.integer(n_per_group),
        length(groups)
      ),
      groups
    )

  } else {

    if (
      is.null(names(n_per_group)) ||
      !all(groups %in% names(n_per_group))
    ) {
      stop(
        "n_per_group must be a single integer or a named vector ",
        "containing low, intermediate, and high.",
        call. = FALSE
      )
    }

    n_per_group <- n_per_group[groups]

    if (
      !is.numeric(n_per_group) ||
      anyNA(n_per_group) ||
      any(!is.finite(n_per_group)) ||
      any(n_per_group < 0) ||
      any(n_per_group != as.integer(n_per_group))
    ) {
      stop(
        "n_per_group must contain non-negative integers.",
        call. = FALSE
      )
    }

    n_per_group <- as.integer(
      n_per_group
    )

    names(n_per_group) <- groups
  }

  if (
    length(seed) != 1L ||
    !is.numeric(seed) ||
    is.na(seed) ||
    !is.finite(seed) ||
    seed != as.integer(seed)
  ) {
    stop(
      "seed must be a single integer.",
      call. = FALSE
    )
  }

  seed <- as.integer(seed)

  if (
    length(strict) != 1L ||
    !is.logical(strict) ||
    is.na(strict)
  ) {
    stop(
      "strict must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  methods <- c(
    low = low_method,
    intermediate = intermediate_method,
    high = high_method
  )

  allowed_methods <- c(
    "ranked",
    "random"
  )

  if (
    anyNA(methods) ||
    !all(methods %in% allowed_methods)
  ) {
    stop(
      "Selection methods must be either 'ranked' or 'random'.",
      call. = FALSE
    )
  }

  candidates <- classified |>
    dplyr::filter(
      as.character(.data$variability_class) %in% groups
    ) |>
    dplyr::mutate(
      variability_class = as.character(
        .data$variability_class
      )
    )

  group_seeds <- c(
    low = seed,
    intermediate = seed + 1L,
    high = seed + 2L
  )

  select_one_group <- function(group_name) {

    x <- candidates |>
      dplyr::filter(
        .data$variability_class == group_name
      )

    n_requested <- n_per_group[[group_name]]
    n_available <- nrow(x)

    if (n_available < n_requested) {

      msg <- paste0(
        "Requested ",
        n_requested,
        " ",
        group_name,
        " taxa, but only ",
        n_available,
        " satisfy the criteria."
      )

      if (strict) {
        stop(
          msg,
          call. = FALSE
        )
      }

      warning(
        msg,
        " Returning all available taxa.",
        call. = FALSE
      )

      n_requested <- n_available
    }

    method <- methods[[group_name]]

    if (n_requested == 0L) {

      x <- x[0, , drop = FALSE]

      x$selection_method <- character()
      x$selection_rank <- integer()

      return(x)
    }

    # Ranked selection --------------------------------------------------

    if (method == "ranked") {

      if (group_name == "low") {

        x <- x |>
          dplyr::arrange(
            dplyr::desc(.data$mean),
            .data$cv,
            dplyr::desc(.data$report_rate),
            .data$taxid
          )

      } else if (group_name == "intermediate") {

        x <- x |>
          dplyr::arrange(
            dplyr::desc(.data$mean),
            abs(.data$cv - 75),
            .data$taxid
          )

      } else if (group_name == "high") {

        x <- x |>
          dplyr::arrange(
            dplyr::desc(.data$mean),
            .data$report_rate,
            dplyr::desc(.data$cv),
            .data$taxid
          )
      }

      x <- dplyr::slice_head(
        x,
        n = n_requested
      )
    }

    # Random selection --------------------------------------------------

    if (method == "random") {

      x <- x |>
        dplyr::arrange(
          .data$taxid
        )

      set.seed(
        group_seeds[[group_name]]
      )

      idx <- sample.int(
        n = nrow(x),
        size = n_requested,
        replace = FALSE
      )

      x <- x[idx, , drop = FALSE]
    }

    x |>
      dplyr::mutate(
        selection_method = method
      )
  }

  selected <- purrr::map_dfr(
    groups,
    select_one_group
  )

  selected <- selected |>
    dplyr::mutate(
      variability_class = factor(
        .data$variability_class,
        levels = groups
      )
    ) |>
    dplyr::group_by(
      .data$variability_class
    ) |>
    dplyr::arrange(
      dplyr::desc(.data$mean),
      .by_group = TRUE
    ) |>
    dplyr::mutate(
      selection_rank = dplyr::row_number()
    ) |>
    dplyr::ungroup()

  selected
}


#' Classify and select taxa for in silico evaluation
#'
#' Convenience wrapper around [classify_inspect_taxa()] and
#' [select_inspect_taxa()].
#'
#' By default, low and high variability taxa are selected by ranking,
#' favoring better represented taxa, while intermediate variability taxa
#' are sampled randomly using a reproducible seed.
#'
#' @inheritParams classify_inspect_taxa
#' @inheritParams select_inspect_taxa
#'
#' @return A list containing:
#'
#' * `classified`: all taxa at the requested taxonomic level with assigned
#'   variability classes.
#' * `candidates`: taxa satisfying low, intermediate, or high criteria.
#' * `selected`: taxa chosen for in silico evaluation.
#' * `summary`: number of candidate and selected taxa in each class.
#'
#' @examples
#' \dontrun{
#' selection <- make_insilico_taxa(
#'   stats,
#'   tax_level = "S",
#'   n_per_group = 10,
#'   seed = 20260917
#' )
#'
#' selection$selected
#' selection$summary
#' }
#'
#' @export
make_insilico_taxa <- function(
    stats,
    tax_level = "S",
    n_per_group = 10,
    seed = 20260917L,

    low_report_min = 0.9,
    low_cv_max = 40,

    mid_report_min = 0.6,
    mid_report_max = 0.8,
    mid_cv_min = 50,
    mid_cv_max = 100,

    high_report_max = 0.3,
    high_cv_min = 160,

    min_mean = 0,

    low_method = "ranked",
    intermediate_method = "random",
    high_method = "ranked",

    strict = FALSE
) {

  classified <- classify_inspect_taxa(
    stats = stats,
    tax_level = tax_level,

    low_report_min = low_report_min,
    low_cv_max = low_cv_max,

    mid_report_min = mid_report_min,
    mid_report_max = mid_report_max,
    mid_cv_min = mid_cv_min,
    mid_cv_max = mid_cv_max,

    high_report_max = high_report_max,
    high_cv_min = high_cv_min,

    min_mean = min_mean
  )

  candidates <- classified |>
    dplyr::filter(
      as.character(.data$variability_class) != "other"
    )

  selected <- select_inspect_taxa(
    classified = classified,
    n_per_group = n_per_group,
    seed = seed,

    low_method = low_method,
    intermediate_method = intermediate_method,
    high_method = high_method,

    strict = strict
  )

  summary <- classified |>
    dplyr::count(
      .data$variability_class,
      name = "n_candidates",
      .drop = FALSE
    ) |>
    dplyr::mutate(
      n_selected = vapply(
        as.character(.data$variability_class),
        function(group) {
          sum(
            as.character(
              selected$variability_class
            ) == group
          )
        },
        integer(1)
      )
    )

  list(
    classified = classified,
    candidates = candidates,
    selected = selected,
    summary = summary
  )
}
