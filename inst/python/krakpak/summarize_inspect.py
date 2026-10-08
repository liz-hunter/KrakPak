#' Summarize Kraken classifications across replicate databases
#'
#' Reads Kraken2 classification outputs for one in silico sample across
#' replicate databases, summarizes classification outcomes, and constructs
#' a per-read classification table.
#'
#' Expected Kraken output filenames follow the pattern
#' `{accession}_{dataset}_db{rep}.txt`.
#'
#' Kraken outputs may contain either numeric taxids or taxon strings
#' produced using `--use-names`, such as
#' `"Olea europaea (taxid 4146)"`.
#'
#' Unclassified reads are retained explicitly as `"unclassified"` in the
#' per-read classification table.
#'
#' @param accession NCBI assembly accession for the evaluated genome.
#' @param dataset Dataset label used in Kraken output filenames.
#' @param true_taxid Expected taxonomic ID for the genome.
#' @param kraken_dir Directory containing Kraken2 output files.
#' @param n_db_reps Number of replicate Kraken databases. Defaults to 10.
#' @param out_dir Optional directory in which to write the classification
#'   and summary tables. If `NULL`, results are returned without writing.
#'
#' @return A list containing:
#'
#' * `classification`: per-read predicted taxids across database replicates.
#' * `summary`: classification metrics for each database replicate.
#'
#' @export
summarize_classification <- function(
    accession,
    dataset,
    true_taxid,
    kraken_dir,
    n_db_reps = 10,
    out_dir = NULL
) {

  # -----------------------------------------------------------------------
  # Input validation
  # -----------------------------------------------------------------------

  if (
    length(accession) != 1L ||
    !is.character(accession) ||
    is.na(accession) ||
    !nzchar(accession)
  ) {
    stop(
      "accession must be a single non-empty character value.",
      call. = FALSE
    )
  }

  if (
    length(dataset) != 1L ||
    !is.character(dataset) ||
    is.na(dataset) ||
    !nzchar(dataset)
  ) {
    stop(
      "dataset must be a single non-empty character value.",
      call. = FALSE
    )
  }

  if (
    length(true_taxid) != 1L ||
    is.na(true_taxid)
  ) {
    stop(
      "true_taxid must contain exactly one taxid.",
      call. = FALSE
    )
  }

  true_taxid <- as.character(true_taxid)

  if (
    length(kraken_dir) != 1L ||
    !is.character(kraken_dir) ||
    is.na(kraken_dir) ||
    !nzchar(kraken_dir) ||
    !dir.exists(kraken_dir)
  ) {
    stop(
      "kraken_dir must be an existing directory.",
      call. = FALSE
    )
  }

  if (
    length(n_db_reps) != 1L ||
    !is.numeric(n_db_reps) ||
    is.na(n_db_reps) ||
    !is.finite(n_db_reps) ||
    n_db_reps < 1 ||
    n_db_reps != as.integer(n_db_reps)
  ) {
    stop(
      "n_db_reps must be a positive integer.",
      call. = FALSE
    )
  }

  n_db_reps <- as.integer(n_db_reps)

  if (!is.null(out_dir)) {

    if (
      length(out_dir) != 1L ||
      !is.character(out_dir) ||
      is.na(out_dir) ||
      !nzchar(out_dir)
    ) {
      stop(
        "out_dir must be a single directory path or NULL.",
        call. = FALSE
      )
    }
  }

  # -----------------------------------------------------------------------
  # Containers
  # -----------------------------------------------------------------------

  classification <- NULL

  summaries <- vector(
    "list",
    n_db_reps
  )

  # -----------------------------------------------------------------------
  # Process database replicates
  # -----------------------------------------------------------------------

  for (db_rep in seq_len(n_db_reps)) {

    filename <- paste0(
      accession,
      "_",
      dataset,
      "_db",
      db_rep,
      ".txt"
    )

    path <- file.path(
      kraken_dir,
      filename
    )

    # ---------------------------------------------------------------------
    # Missing Kraken file
    # ---------------------------------------------------------------------

    if (!file.exists(path)) {

      summaries[[db_rep]] <- tibble::tibble(
        accession = accession,
        dataset = dataset,
        db_rep = db_rep,
        true_taxid = true_taxid,

        n_reads = NA_integer_,
        n_classified = NA_integer_,
        n_correct = NA_integer_,
        n_misclassified = NA_integer_,
        n_unclassified = NA_integer_,

        pct_correct_all = NA_real_,
        pct_misclassified_all = NA_real_,
        pct_unclassified = NA_real_,
        pct_correct_classified = NA_real_,

        file_found = FALSE
      )

      next
    }

    # ---------------------------------------------------------------------
    # Read Kraken output
    # ---------------------------------------------------------------------

    x <- read_kraken_classification(
      path
    )

    # ---------------------------------------------------------------------
    # Classification outcomes
    # ---------------------------------------------------------------------

    n_reads <- nrow(x)

    classified <- x$status == "C"

    unclassified <- x$status == "U"

    correct <- classified &
      x$pred_taxid == true_taxid

    misclassified <- classified &
      x$pred_taxid != true_taxid

    n_classified <- sum(
      classified
    )

    n_correct <- sum(
      correct
    )

    n_misclassified <- sum(
      misclassified
    )

    n_unclassified <- sum(
      unclassified
    )

    # Every read should belong to exactly one outcome category.
    if (
      n_correct +
        n_misclassified +
        n_unclassified != n_reads
    ) {
      stop(
        "Classification outcome counts do not sum to total reads in ",
        path,
        ".",
        call. = FALSE
      )
    }

    # ---------------------------------------------------------------------
    # Percentages
    # ---------------------------------------------------------------------

    pct_correct_all <- if (n_reads > 0L) {
      100 * n_correct / n_reads
    } else {
      NA_real_
    }

    pct_misclassified_all <- if (n_reads > 0L) {
      100 * n_misclassified / n_reads
    } else {
      NA_real_
    }

    pct_unclassified <- if (n_reads > 0L) {
      100 * n_unclassified / n_reads
    } else {
      NA_real_
    }

    pct_correct_classified <- if (n_classified > 0L) {
      100 * n_correct / n_classified
    } else {
      NA_real_
    }

    # ---------------------------------------------------------------------
    # Summary row
    # ---------------------------------------------------------------------

    summaries[[db_rep]] <- tibble::tibble(
      accession = accession,
      dataset = dataset,
      db_rep = db_rep,
      true_taxid = true_taxid,

      n_reads = n_reads,
      n_classified = n_classified,
      n_correct = n_correct,
      n_misclassified = n_misclassified,
      n_unclassified = n_unclassified,

      pct_correct_all = pct_correct_all,
      pct_misclassified_all = pct_misclassified_all,
      pct_unclassified = pct_unclassified,
      pct_correct_classified = pct_correct_classified,

      file_found = TRUE
    )

    # ---------------------------------------------------------------------
    # Per-read classifications
    # ---------------------------------------------------------------------

    column_name <- paste0(
      dataset,
      "_db",
      db_rep
    )

    current <- x |>
      dplyr::select(
        .data$read_id,
        .data$pred_taxid
      )

    names(current)[2] <- column_name

    if (is.null(classification)) {

      classification <- current

    } else {

      classification <- dplyr::full_join(
        classification,
        current,
        by = "read_id",
        relationship = "one-to-one"
      )
    }
  }

  # -----------------------------------------------------------------------
  # Require at least one Kraken result
  # -----------------------------------------------------------------------

  if (is.null(classification)) {
    stop(
      "No Kraken output files were found for ",
      accession,
      " (",
      dataset,
      ").",
      call. = FALSE
    )
  }

  # -----------------------------------------------------------------------
  # Add true taxid
  # -----------------------------------------------------------------------

  classification <- classification |>
    dplyr::mutate(
      true_taxid = true_taxid,
      .after = .data$read_id
    )

  summary <- dplyr::bind_rows(
    summaries
  )

  result <- list(
    classification = classification,
    summary = summary
  )

  # -----------------------------------------------------------------------
  # Optional output
  # -----------------------------------------------------------------------

  if (!is.null(out_dir)) {

    dir.create(
      out_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    classification_file <- file.path(
      out_dir,
      paste0(
        accession,
        "_",
        dataset,
        "_classification.tsv"
      )
    )

    summary_file <- file.path(
      out_dir,
      paste0(
        accession,
        "_",
        dataset,
        "_summary.tsv"
      )
    )

    readr::write_tsv(
      classification,
      classification_file
    )

    readr::write_tsv(
      summary,
      summary_file
    )

    result$classification_file <- classification_file
    result$summary_file <- summary_file
  }

  result
}


# -------------------------------------------------------------------------
# Internal helper: extract numeric taxid from Kraken taxon field
# -------------------------------------------------------------------------

.read_kraken_taxid <- function(x) {

  x <- trimws(
    as.character(x)
  )

  out <- rep(
    NA_character_,
    length(x)
  )

  # Already numeric
  numeric_taxid <- grepl(
    "^[0-9]+$",
    x
  )

  out[numeric_taxid] <- x[numeric_taxid]

  # Kraken --use-names format:
  # "Olea europaea (taxid 4146)"
  named_taxid <- stringr::str_match(
    x,
    "\\(taxid\\s+([0-9]+)\\)\\s*$"
  )

  has_named_taxid <- !is.na(
    named_taxid[, 2]
  )

  out[has_named_taxid] <- named_taxid[
    has_named_taxid,
    2
  ]

  out
}


# -------------------------------------------------------------------------
# Internal helper: read one Kraken classification file
# -------------------------------------------------------------------------

read_kraken_classification <- function(path) {

  if (!file.exists(path)) {
    stop(
      "Kraken classification file does not exist: ",
      path,
      call. = FALSE
    )
  }

  x <- readr::read_tsv(
    path,
    col_names = FALSE,
    col_select = 1:3,
    col_types = "ccc",
    show_col_types = FALSE
  )

  names(x) <- c(
    "status",
    "read_id",
    "pred_taxid"
  )

  parsed_taxid <- .read_kraken_taxid(
    x$pred_taxid
  )

  # Classified reads should always contain a recognizable taxid.
  bad <- x$status == "C" &
    is.na(parsed_taxid)

  if (any(bad)) {
    stop(
      sum(bad),
      " classified read(s) in ",
      path,
      " did not contain a recognizable taxid.",
      call. = FALSE
    )
  }

  # Explicitly distinguish unclassified reads from missing values.
  x$pred_taxid <- dplyr::case_when(
    x$status == "U" ~ "unclassified",
    x$status == "C" ~ parsed_taxid,
    TRUE ~ NA_character_
  )

  x
}
