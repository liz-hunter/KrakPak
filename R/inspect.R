#' Read and filter Kraken2 inspect reports
#'
#' Reads one or more Kraken2 inspect reports from a directory, combines
#' them into a single tibble, and adds a database identifier derived from
#' each filename. Exact minimizer percentages are recalculated from the
#' inclusive minimizer counts because Kraken2 reports rounded percentages.
#'
#' Files are selected by matching filenames that begin with `prefix`.
#'
#' @param directory Directory containing Kraken2 inspect reports.
#' @param prefix Character string identifying files to read. Files must
#'   begin with this string.
#' @param level Taxonomic level to retain. One of `"species"`, `"genus"`,
#'   `"family"`, `"order"`, `"class"`, `"phylum"`, `"kingdom"`, `"domain"`,
#'   `"unclassified"`, or `"all"`. Defaults to `"species"`.
#' @param set_taxids Optional vector of taxids, or data frame containing a
#'   column named `taxid`. If supplied, only these taxids are retained.
#' @param exclude_taxids Optional vector of taxids, or data frame containing
#'   a column named `taxid`. If supplied, these taxids are removed.
#'   A taxid may not appear in both `set_taxids` and `exclude_taxids`.
#' @param fuzzy Logical. Controls whether Kraken2 numbered intermediate ranks
#'   (e.g. `G1`, `G2`, `S1`, `S2`) are retained. If `FALSE`, only canonical
#'   ranks are retained. If `TRUE`, intermediate ranks are also retained.
#'   When `level = "all"` and `fuzzy = TRUE`, all taxonomic ranks are
#'   returned. Defaults to `FALSE`.
#'
#' @return A tibble containing the combined and filtered Kraken2 inspect
#'   reports. The original rounded `perc_comp` column is retained and
#'   `perc_comp_exact` contains the recalculated percentage.
#'
#' @export
#'
read_inspect <- function(
    directory,
    prefix,
    level = "species",
    fuzzy = FALSE,
    set_taxids = NULL,
    exclude_taxids = NULL
) {

  # Column structure of a Kraken2 inspect report
  report_head <- c(
    "perc_comp",
    "incl_min_count",
    "excl_min_count",
    "level",
    "taxid",
    "name"
  )

  # Translate human-readable taxonomic ranks to Kraken rank codes
  rank_codes <- c(
    unclassified = "U",
    root         = "R",
    domain       = "D",
    kingdom      = "K",
    phylum       = "P",
    class        = "C",
    order        = "O",
    family       = "F",
    genus        = "G",
    species      = "S"
  )

  # Check directory
  if (!dir.exists(directory)) {
    stop("Directory does not exist: ", directory)
  }

  # Check level
  level <- tolower(level)

  if (!level %in% c(names(rank_codes), "all")) {
    stop(
      "`level` must be one of: ",
      paste(c(names(rank_codes), "all"), collapse = ", ")
    )
  }

  # Find files beginning with prefix
  files <- list.files(
    path = directory,
    full.names = TRUE
  )

  files <- files[
    startsWith(basename(files), prefix)
  ]

  if (length(files) == 0) {
    stop(
      "No files beginning with '",
      prefix,
      "' found in ",
      directory
    )
  }

  # Read individual inspect reports
  inspect <- lapply(files, function(file) {

    dat <- readr::read_tsv(
      file,
      col_names = report_head,
      show_col_types = FALSE,
      trim_ws = FALSE
    )

    if (nrow(dat) == 0) {
      stop("Inspect report is empty: ", file)
    }

    # Preserve Kraken report order
    dat$report_row <- seq_len(nrow(dat))

    # Record taxonomic indentation without modifying the original name
    dat$indent <- nchar(dat$name) -
      nchar(sub("^ +", "", dat$name))

    # Database identifier from filename
    dat$db <- sub(
      "_inspect\\.txt$",
      "",
      basename(file)
    )

    # Root inclusive count represents all minimizers in the database
    root_count <- dat$incl_min_count[
      dat$level == "R" & dat$taxid == 1
    ]

    if (length(root_count) != 1) {
      stop(
        "Could not identify exactly one root row in: ",
        file
      )
    }

    # Recalculate percentage without Kraken's rounding
    dat$perc_comp_exact <-
      (dat$incl_min_count / root_count) * 100

    dat
  })

  # Combine reports into one tibble
  inspect <- dplyr::bind_rows(inspect)

  # Normalize set_taxids to a vector
  if (!is.null(set_taxids)) {

    if (is.data.frame(set_taxids)) {

      if (!"taxid" %in% names(set_taxids)) {
        stop(
          "`set_taxids` data frame must contain a column named `taxid`."
        )
      }

      set_taxids <- set_taxids$taxid
    }

    set_taxids <- unique(set_taxids)
  }

  # Normalize exclude_taxids to a vector
  if (!is.null(exclude_taxids)) {

    if (is.data.frame(exclude_taxids)) {

      if (!"taxid" %in% names(exclude_taxids)) {
        stop(
          "`exclude_taxids` data frame must contain a column named `taxid`."
        )
      }

      exclude_taxids <- exclude_taxids$taxid
    }

    exclude_taxids <- unique(exclude_taxids)
  }

  # Fail if any taxid is both included and excluded
  if (!is.null(set_taxids) && !is.null(exclude_taxids)) {

    overlap <- intersect(set_taxids, exclude_taxids)

    if (length(overlap) > 0) {
      stop(
        "The following taxids are present in both `set_taxids` and ",
        "`exclude_taxids`: ",
        paste(overlap, collapse = ", ")
      )
    }
  }

  # Optionally restrict to specified taxids
  if (!is.null(set_taxids)) {
    inspect <- inspect[
      inspect$taxid %in% set_taxids,
    ]
  }

  # Optionally remove specified taxids
  if (!is.null(exclude_taxids)) {
    inspect <- inspect[
      !inspect$taxid %in% exclude_taxids,
    ]
  }

  # Filter by taxonomic level and fuzzy rank behavior
  if (level == "all") {

    if (!fuzzy) {

      # Keep canonical ranks only
      canonical_codes <- c("U", "R", unname(rank_codes))

      inspect <- inspect[
        inspect$level %in% canonical_codes,
      ]
    }

    # If level == "all" and fuzzy == TRUE,
    # retain every rank in the report

  } else {

    rank_code <- unname(rank_codes[[level]])

    if (fuzzy) {

      # Canonical rank plus numbered intermediate ranks
      # e.g. S, S1, S2, ...
      keep <- grepl(
        paste0("^", rank_code, "[0-9]*$"),
        inspect$level
      )

    } else {

      # Exact canonical rank only
      keep <- inspect$level == rank_code
    }

    inspect <- inspect[keep, ]
  }

inspect
}
