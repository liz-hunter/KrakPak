#' Extract descendants of a taxon from a Kraken2 inspect report
#'
#' Extracts the descendants of a specified taxid from a Kraken2 inspect
#' report using the row order and indentation encoded in the report.
#'
#' The input should contain the complete Kraken taxonomy tree, typically
#' produced with `read_inspect(level = "all", fuzzy = TRUE)` and without
#' taxid inclusion or exclusion filters.
#'
#' @param inspect A tibble returned by `read_inspect()`.
#' @param taxid A single parent taxid whose descendants should be extracted.
#' @param include_parent Logical. If `TRUE`, include the parent taxon itself
#'   in the returned tibble. Defaults to `TRUE`.
#' @param immediate Logical. If `TRUE`, return only the immediate children
#'   of the parent. If `FALSE`, return all descendants. Defaults to `FALSE`.
#' @param validate Logical. If `TRUE`, validate the extracted hierarchy
#'   against Kraken2 minimizer counts. Defaults to `TRUE`.
#'
#' @return A tibble containing the requested taxon and/or its descendants
#'   for each database in which the taxid occurs.
#'
#' @export
inspect_children <- function(
    inspect,
    taxid,
    include_parent = TRUE,
    immediate = FALSE,
    validate = TRUE
) {

  required_cols <- c(
    "db",
    "report_row",
    "indent",
    "taxid",
    "name",
    "incl_min_count",
    "excl_min_count"
  )

  missing_cols <- setdiff(required_cols, names(inspect))

  if (length(missing_cols) > 0) {
    stop(
      "`inspect` is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (length(taxid) != 1 || is.na(taxid)) {
    stop("`taxid` must contain exactly one non-missing taxid.")
  }

  if (!is.logical(include_parent) ||
      length(include_parent) != 1 ||
      is.na(include_parent)) {
    stop("`include_parent` must be either TRUE or FALSE.")
  }

  if (!is.logical(immediate) ||
      length(immediate) != 1 ||
      is.na(immediate)) {
    stop("`immediate` must be either TRUE or FALSE.")
  }

  if (!is.logical(validate) ||
      length(validate) != 1 ||
      is.na(validate)) {
    stop("`validate` must be either TRUE or FALSE.")
  }

  # Process each database independently
  dbs <- split(inspect, inspect$db)

  out <- lapply(dbs, function(dat) {

    dat <- dat[order(dat$report_row), ]

    parent_idx <- which(dat$taxid == taxid)

    # Taxon may not occur in every database
    if (length(parent_idx) == 0) {
      return(NULL)
    }

    if (length(parent_idx) > 1) {
      stop(
        "Taxid ",
        taxid,
        " occurs more than once in database `",
        dat$db[1],
        "`."
      )
    }

    parent_idx <- parent_idx[1]
    parent_indent <- dat$indent[parent_idx]

    # Find the end of this taxon's subtree.
    # Descendants continue until indentation returns to the
    # parent's depth or higher.
    if (parent_idx == nrow(dat)) {

      descendant_idx <- integer(0)

    } else {

      following_idx <- seq.int(parent_idx + 1, nrow(dat))

      subtree_end <- which(
        dat$indent[following_idx] <= parent_indent
      )

      if (length(subtree_end) == 0) {
        last_idx <- nrow(dat)
      } else {
        last_idx <- following_idx[subtree_end[1]] - 1
      }

      if (last_idx >= parent_idx + 1) {
        descendant_idx <- seq.int(parent_idx + 1, last_idx)
      } else {
        descendant_idx <- integer(0)
      }
    }

    # Restrict to immediate children if requested
    if (immediate && length(descendant_idx) > 0) {

      child_indent <- min(dat$indent[descendant_idx])

      descendant_idx <- descendant_idx[
        dat$indent[descendant_idx] == child_indent
      ]
    }

    # Validate tree reconstruction against minimizer counts
    if (validate) {

      expected <- dat$incl_min_count[parent_idx]

      if (immediate) {

        # Parent-exclusive minimizers plus inclusive minimizers
        # from each immediate child should reconstruct the parent.
        observed <-
          dat$excl_min_count[parent_idx] +
          sum(dat$incl_min_count[descendant_idx])

      } else {

        # Exclusive counts across the entire subtree should
        # reconstruct the parent's inclusive count.
        subtree_idx <- c(parent_idx, descendant_idx)

        observed <- sum(dat$excl_min_count[subtree_idx])
      }

      if (!isTRUE(all.equal(
        observed,
        expected,
        tolerance = 0
      ))) {
        stop(
          "Tree validation failed for taxid ",
          taxid,
          " in database `",
          dat$db[1],
          "`. Expected ",
          expected,
          " minimizers but reconstructed ",
          observed,
          ". Make sure the input contains the complete hierarchy ",
          "(e.g. `level = \"all\"`, `fuzzy = TRUE`, and no taxid filters)."
        )
      }
    }

    keep_idx <- descendant_idx

    if (include_parent) {
      keep_idx <- c(parent_idx, keep_idx)
    }

    dat[keep_idx, , drop = FALSE]
  })

  out <- out[!vapply(out, is.null, logical(1))]

  if (length(out) == 0) {
    stop(
      "Taxid ",
      taxid,
      " was not found in any database."
    )
  }

  dplyr::bind_rows(out)
}
