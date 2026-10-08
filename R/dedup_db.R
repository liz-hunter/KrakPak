# db_prep.R
#
# Functions for preparing NCBI genome metadata for
# custom Kraken2 database construction


# Internal constants -----------------------------------------------------

.ncbi_required_columns <- c(
  "assembly_name",
  "assembly_accession",
  "assembly_paired_assembly_accession",
  "organism_name",
  "organism_taxonomic_id"
)

.valid_assembly_levels <- c(
  "Complete Genome",
  "Chromosome",
  "Scaffold",
  "Contig"
)

# Internal helpers -------------------------------------------------------

.validate_threshold <- function(x, arg, max = Inf) {
  if (is.null(x)) {
    return(NULL)
  }

  if (
    length(x) != 1L ||
    is.na(x) ||
    !is.numeric(x) ||
    !is.finite(x) ||
    x < 0 ||
    x > max
  ) {
    if (is.finite(max)) {
      stop(
        arg,
        " must be a single numeric value between 0 and ",
        max,
        ".",
        call. = FALSE
      )
    }

    stop(
      arg,
      " must be a single non-negative numeric value.",
      call. = FALSE
    )
  }

  x
}


.require_column <- function(df, column, filter_name) {
  if (!column %in% names(df)) {
    stop(
      "Cannot filter by ",
      filter_name,
      ": column '",
      column,
      "' is not present in the input.",
      call. = FALSE
    )
  }
}


.as_numeric_quality <- function(x) {
  if (is.numeric(x)) {
    return(x)
  }

  suppressWarnings(
    readr::parse_number(
      as.character(x),
      na = c("", "NA", "N/A", "na", "n/a", "-")
    )
  )
}


.normalize_level <- function(x) {
  x <- stringr::str_to_lower(
    stringr::str_trim(as.character(x))
  )

  # Allow the convenient shorthand "complete".
  x[x == "complete"] <- "complete genome"

  x
}


#' Remove redundant GenBank/RefSeq assemblies
#'
#' Reads an NCBI genome metadata table and removes redundant paired
#' GenBank/RefSeq assemblies. When both versions of an assembly are
#' present, the RefSeq (`GCF_`) assembly is retained.
#'
#' The input must contain `assembly_name`, `assembly_accession`,
#' `assembly_paired_assembly_accession`, `organism_name`, and
#' `organism_taxonomic_id`. Additional metadata columns are retained.
#'
#' Standardized columns named `accession`, `name`, `taxid`, and `filename`
#' are added for downstream KrakPak functions.
#'
#' @param infile Path to an NCBI genome metadata TSV.
#' @param report Logical; print a summary of duplicate removal.
#'
#' @return A tibble containing one row per retained assembly, with the
#'   original NCBI metadata plus standardized KrakPak columns.
#'
#' @examples
#' \dontrun{
#' db <- dedup_db("ncbi_genomes.tsv")
#' }
#'
#' @export
dedup_db <- function(infile, report = TRUE) {

  if (
    !is.character(infile) ||
    length(infile) != 1L ||
    is.na(infile) ||
    infile == ""
  ) {
    stop(
      "infile must be a path to an NCBI genome metadata TSV.",
      call. = FALSE
    )
  }

  if (!file.exists(infile)) {
    stop(
      "Input file does not exist: ",
      infile,
      call. = FALSE
    )
  }

  db <- readr::read_tsv(
    infile,
    show_col_types = FALSE,
    progress = FALSE
  ) |>
    janitor::clean_names()

  missing_cols <- setdiff(
    .ncbi_required_columns,
    names(db)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Check for required values
  required_values <- c(
    "assembly_name",
    "assembly_accession",
    "organism_name",
    "organism_taxonomic_id"
  )

  for (column in required_values) {
    value <- db[[column]]

    bad <- is.na(value) |
      stringr::str_trim(as.character(value)) == ""

    if (any(bad)) {
      stop(
        "Column '",
        column,
        "' contains missing or empty values (n=",
        sum(bad),
        ").",
        call. = FALSE
      )
    }
  }

  # Normalize empty paired-accession values
  db <- db |>
    dplyr::mutate(
      assembly_paired_assembly_accession =
        stringr::str_trim(
          assembly_paired_assembly_accession
        ),
      assembly_paired_assembly_accession =
        dplyr::na_if(
          assembly_paired_assembly_accession,
          ""
        )
    )

  n_input <- nrow(db)

  # Pair the GCA/GCF accessions that share the same assembly identifier
  db_ranked <- db |>
    dplyr::mutate(
      .key_source = dplyr::coalesce(
        assembly_paired_assembly_accession,
        assembly_accession
      ),
      .pair_key = .key_source |>
        stringr::str_remove("^GC[AF]_") |>
        stringr::str_remove("\\..*$"),
      .is_refseq = stringr::str_starts(
        assembly_accession,
        "GCF_"
      )
    ) |>
    dplyr::group_by(.pair_key) |>
    dplyr::arrange(
      dplyr::desc(.is_refseq),
      assembly_accession,
      .by_group = TRUE
    ) |>
    dplyr::mutate(
      .rank = dplyr::row_number()
    ) |>
    dplyr::ungroup()

  removed <- db_ranked |>
    dplyr::filter(.rank > 1L)

  kept <- db_ranked |>
    dplyr::filter(.rank == 1L) |>
    dplyr::select(
      -.key_source,
      -.pair_key,
      -.is_refseq,
      -.rank
    )|>
    dplyr::mutate(
      accession = assembly_accession,
      name = organism_name,
      taxid = organism_taxonomic_id,
      filename = paste0(
        assembly_accession,
        "_",
        assembly_name,
        "_genomic.fna"
      )
    ) |>
    dplyr::relocate(
      accession,
      name,
      taxid,
      filename
    )

  if (report) {
    cat(
      "\nDatabase de-duplication summary:\n",
      "  Input assemblies:           ", n_input, "\n",
      "  Redundant assemblies:       ", nrow(removed), "\n",
      "  Assemblies retained:        ", nrow(kept), "\n\n",
      sep = ""
    )
  }

  kept
}

#' Filter genome assemblies by assembly quality
#'
#' Applies optional quality filters to an assembly table produced by
#' [dedup_db()]. Filters are applied only when explicitly requested.
#'
#' Missing quality values are retained by default. Set
#' `keep_missing = FALSE` to remove assemblies lacking a value for any
#' requested quality metric.
#'
#' @param db A data frame produced by [dedup_db()].
#' @param assembly_level Optional character vector of NCBI assembly levels
#'   to retain. Valid values are `"Complete Genome"`, `"Chromosome"`,
#'   `"Scaffold"`, and `"Contig"`.
#' @param min_contig_n50 Optional minimum Contig N50 in base pairs.
#' @param min_scaffold_n50 Optional minimum Scaffold N50 in base pairs.
#' @param min_checkm_completeness Optional minimum CheckM completeness
#'   percentage.
#' @param max_checkm_contamination Optional maximum CheckM contamination
#'   percentage.
#' @param keep_missing Logical; retain assemblies with missing values for
#'   requested quality metrics. Defaults to `TRUE`.
#' @param report Logical; print a summary of filtering.
#'
#' @return A filtered tibble.
#'
#' @examples
#' \dontrun{
#' db <- dedup_db("ncbi_genomes.tsv")
#'
#' complete <- filter_db(
#'   db,
#'   assembly_level = "Complete Genome"
#' )
#'
#' high_quality <- filter_db(
#'   db,
#'   assembly_level = c("Complete Genome", "Chromosome"),
#'   min_contig_n50 = 50000
#' )
#'
#' prokaryotes <- filter_db(
#'   db,
#'   min_checkm_completeness = 90,
#'   max_checkm_contamination = 5
#' )
#' }
#'
#' @export
filter_db <- function(
    db,
    assembly_level = NULL,
    min_contig_n50 = NULL,
    min_scaffold_n50 = NULL,
    min_checkm_completeness = NULL,
    max_checkm_contamination = NULL,
    keep_missing = TRUE,
    report = TRUE
) {

  if (!is.data.frame(db)) {
    stop(
      "db must be a data frame produced by dedup_db().",
      call. = FALSE
    )
  }

  needed <- c(
    "accession",
    "name",
    "taxid",
    "filename"
  )

  missing_cols <- setdiff(needed, names(db))

  if (length(missing_cols) > 0L) {
    stop(
      "db is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      "\nRun dedup_db() before filter_db().",
      call. = FALSE
    )
  }

  if (
    length(keep_missing) != 1L ||
    is.na(keep_missing) ||
    !is.logical(keep_missing)
  ) {
    stop(
      "keep_missing must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  min_contig_n50 <- .validate_threshold(
    min_contig_n50,
    "min_contig_n50"
  )

  min_scaffold_n50 <- .validate_threshold(
    min_scaffold_n50,
    "min_scaffold_n50"
  )

  min_checkm_completeness <- .validate_threshold(
    min_checkm_completeness,
    "min_checkm_completeness",
    max = 100
  )

  max_checkm_contamination <- .validate_threshold(
    max_checkm_contamination,
    "max_checkm_contamination",
    max = 100
  )

  out <- tibble::as_tibble(db)

  n_input <- nrow(out)
  filter_report <- character()


  # Assembly level -------------------------------------------------------

  if (!is.null(assembly_level)) {

    .require_column(
      out,
      "assembly_level",
      "assembly level"
    )

    if (
      !is.character(assembly_level) ||
      length(assembly_level) == 0L ||
      anyNA(assembly_level)
    ) {
      stop(
        "assembly_level must contain one or more assembly levels.",
        call. = FALSE
      )
    }

    assembly_level <- unique(assembly_level)

    bad_levels <- setdiff(
      assembly_level,
      .valid_assembly_levels
    )

    if (length(bad_levels) > 0L) {
      stop(
        "Unknown assembly level(s): ",
        paste(bad_levels, collapse = ", "),
        "\nAllowed: ",
        paste(.valid_assembly_levels, collapse = ", "),
        call. = FALSE
      )
    }

    keep <- out$assembly_level %in% assembly_level
    keep[is.na(out$assembly_level)] <- keep_missing

    n_before <- nrow(out)
    out <- out[keep, , drop = FALSE]

    filter_report <- c(
      filter_report,
      paste0(
        "Assembly level (",
        paste(assembly_level, collapse = ", "),
        "): ",
        n_before - nrow(out),
        " removed"
      )
    )
  }


  # Contig N50 -----------------------------------------------------------

  if (!is.null(min_contig_n50)) {

    .require_column(
      out,
      "contig_n50",
      "Contig N50"
    )

    value <- .as_numeric_quality(
      out$contig_n50
    )

    out$contig_n50 <- value

    keep <- value >= min_contig_n50
    keep[is.na(keep)] <- keep_missing

    n_before <- nrow(out)
    out <- out[keep, , drop = FALSE]

    filter_report <- c(
      filter_report,
      paste0(
        "Contig N50 >= ",
        format(min_contig_n50, big.mark = ","),
        ": ",
        n_before - nrow(out),
        " removed"
      )
    )
  }


  # Scaffold N50 ---------------------------------------------------------

  if (!is.null(min_scaffold_n50)) {

    .require_column(
      out,
      "scaffold_n50",
      "Scaffold N50"
    )

    value <- .as_numeric_quality(
      out$scaffold_n50
    )

    out$scaffold_n50 <- value

    keep <- value >= min_scaffold_n50
    keep[is.na(keep)] <- keep_missing

    n_before <- nrow(out)
    out <- out[keep, , drop = FALSE]

    filter_report <- c(
      filter_report,
      paste0(
        "Scaffold N50 >= ",
        format(min_scaffold_n50, big.mark = ","),
        ": ",
        n_before - nrow(out),
        " removed"
      )
    )
  }


  # CheckM completeness --------------------------------------------------

  if (!is.null(min_checkm_completeness)) {

    .require_column(
      out,
      "check_m_completeness",
      "CheckM completeness"
    )

    value <- .as_numeric_quality(
      out$check_m_completeness
    )

    out$check_m_completeness <- value

    keep <- value >= min_checkm_completeness
    keep[is.na(keep)] <- keep_missing

    n_before <- nrow(out)
    out <- out[keep, , drop = FALSE]

    filter_report <- c(
      filter_report,
      paste0(
        "CheckM completeness >= ",
        min_checkm_completeness,
        "%: ",
        n_before - nrow(out),
        " removed"
      )
    )
  }


  # CheckM contamination -------------------------------------------------

  if (!is.null(max_checkm_contamination)) {

    .require_column(
      out,
      "check_m_contamination",
      "CheckM contamination"
    )

    value <- .as_numeric_quality(
      out$check_m_contamination
    )

    out$check_m_contamination <- value

    keep <- value <= max_checkm_contamination
    keep[is.na(keep)] <- keep_missing

    n_before <- nrow(out)
    out <- out[keep, , drop = FALSE]

    filter_report <- c(
      filter_report,
      paste0(
        "CheckM contamination <= ",
        max_checkm_contamination,
        "%: ",
        n_before - nrow(out),
        " removed"
      )
    )
  }


  # Report ---------------------------------------------------------------

  if (report) {

    cat(
      "\nDatabase quality filtering summary:\n",
      "  Input assemblies:  ", n_input, "\n",
      sep = ""
    )

    if (length(filter_report) == 0L) {
      cat("  No quality filters requested.\n")
    } else {
      for (x in filter_report) {
        cat("  ", x, "\n", sep = "")
      }
    }

    cat(
      "  Missing values:    ",
      if (keep_missing) "retained" else "removed",
      "\n",
      "  Output assemblies: ",
      nrow(out),
      "\n\n",
      sep = ""
    )
  }

  out
}
