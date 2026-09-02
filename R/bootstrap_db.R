# bootstrap.R
#
# Functions for generating bootstrap samples of genome assemblies.


# Internal constants -----------------------------------------------------

.bootstrap_outputs <- c(
  "composite",
  "files",
  "unique_files",
  "accession_counts",
  "taxid_counts",
  "all_accessions"
)


# Internal helpers -------------------------------------------------------

.sanitize_bootstrap_label <- function(label) {

  if (
    !is.character(label) ||
    length(label) != 1L ||
    is.na(label) ||
    stringr::str_trim(label) == ""
  ) {
    stop(
      "label must be a single non-empty character value.",
      call. = FALSE
    )
  }

  label_safe <- label |>
    stringr::str_trim() |>
    stringr::str_replace_all("\\s+", "_") |>
    stringr::str_replace_all("[^A-Za-z0-9_-]", "_")

  if (label_safe == "") {
    stop(
      "label contains no usable characters.",
      call. = FALSE
    )
  }

  label_safe
}


.validate_bootstrap_outputs <- function(outputs) {

  if (
    !is.character(outputs) ||
    length(outputs) == 0L ||
    anyNA(outputs)
  ) {
    stop(
      "outputs must contain one or more output types.",
      call. = FALSE
    )
  }

  bad <- setdiff(
    outputs,
    .bootstrap_outputs
  )

  if (length(bad) > 0L) {
    stop(
      "Unknown output(s): ",
      paste(bad, collapse = ", "),
      "\nAllowed: ",
      paste(.bootstrap_outputs, collapse = ", "),
      call. = FALSE
    )
  }

  unique(outputs)
}


.validate_bootstrap_db <- function(db) {

  if (!is.data.frame(db)) {
    stop(
      "db must be a data frame produced by dedup_db() or filter_db().",
      call. = FALSE
    )
  }

  needed <- c(
    "accession",
    "name",
    "taxid",
    "filename"
  )

  missing_cols <- setdiff(
    needed,
    names(db)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "db is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      "\nRun dedup_db() before bootstrap_db().",
      call. = FALSE
    )
  }

  if (nrow(db) == 0L) {
    stop(
      "db contains no assemblies.",
      call. = FALSE
    )
  }

  for (column in needed) {

    value <- db[[column]]

    bad <- is.na(value) |
      stringr::str_trim(as.character(value)) == ""

    if (any(bad)) {
      stop(
        "Column '",
        column,
        "' contains missing or empty values.",
        call. = FALSE
      )
    }
  }

  if (anyDuplicated(db$accession)) {
    duplicates <- unique(
      db$accession[duplicated(db$accession)]
    )

    stop(
      "db contains duplicated accessions: ",
      paste(utils::head(duplicates, 10L), collapse = ", "),
      if (length(duplicates) > 10L) " ..." else "",
      "\nRun dedup_db() before bootstrap_db().",
      call. = FALSE
    )
  }

  tibble::as_tibble(db)
}


.validate_n_boot <- function(n_boot) {

  if (
    length(n_boot) != 1L ||
    is.na(n_boot) ||
    !is.numeric(n_boot) ||
    !is.finite(n_boot) ||
    n_boot <= 0 ||
    n_boot != floor(n_boot)
  ) {
    stop(
      "n_boot must be a positive integer.",
      call. = FALSE
    )
  }

  as.integer(n_boot)
}


.validate_bootstrap_seed <- function(seed) {

  if (is.null(seed)) {
    return(NULL)
  }

  if (
    length(seed) != 1L ||
    is.na(seed) ||
    !is.numeric(seed) ||
    !is.finite(seed) ||
    seed != floor(seed)
  ) {
    stop(
      "seed must be a single integer.",
      call. = FALSE
    )
  }

  as.integer(seed)
}


.make_bootstrap_matrices <- function(
    db,
    n_boot,
    seed = NULL
) {

  n <- nrow(db)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Each column contains one bootstrap sample of n assemblies,
  # sampled with replacement.
  idx <- matrix(
    sample.int(
      n,
      size = n * n_boot,
      replace = TRUE
    ),
    nrow = n,
    ncol = n_boot
  )

  mat_acc <- matrix(
    db$accession[idx],
    nrow = n,
    ncol = n_boot
  )

  mat_tax <- matrix(
    db$taxid[idx],
    nrow = n,
    ncol = n_boot
  )

  list(
    accessions = mat_acc,
    taxids = mat_tax,
    indices = idx
  )
}


.set_bootstrap_names <- function(mat, label) {

  colnames(mat) <- paste0(
    label,
    "_db",
    seq_len(ncol(mat))
  )

  mat
}


.make_bootstrap_outdir <- function(outdir) {

  if (!dir.exists(outdir)) {

    created <- dir.create(
      outdir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    if (!created && !dir.exists(outdir)) {
      stop(
        "Could not create output directory: ",
        outdir,
        call. = FALSE
      )
    }
  }

  normalizePath(
    outdir,
    winslash = "/",
    mustWork = TRUE
  )
}


.write_bootstrap_composites <- function(
    mat_acc,
    mat_tax,
    outdir,
    label
) {

  acc_path <- file.path(
    outdir,
    paste0(
      "composite_",
      label,
      "_samples.tsv"
    )
  )

  tax_path <- file.path(
    outdir,
    paste0(
      "composite_",
      label,
      "_taxids.tsv"
    )
  )

  readr::write_tsv(
    as.data.frame(mat_acc),
    acc_path
  )

  readr::write_tsv(
    as.data.frame(mat_tax),
    tax_path
  )

  list(
    accessions = acc_path,
    taxids = tax_path
  )
}


.write_bootstrap_files <- function(
    mat_acc,
    outdir,
    label
) {

  n_boot <- ncol(mat_acc)
  paths <- character(n_boot)

  for (i in seq_len(n_boot)) {

    path <- file.path(
      outdir,
      paste0(
        label,
        "_db",
        i,
        ".txt"
      )
    )

    # Duplicates are intentionally preserved because this represents
    # the raw bootstrap draw.
    writeLines(
      mat_acc[, i],
      con = path
    )

    paths[i] <- path
  }

  paths
}


.write_unique_bootstrap_files <- function(
    mat_acc,
    outdir,
    label
) {

  n_boot <- ncol(mat_acc)
  paths <- character(n_boot)

  for (i in seq_len(n_boot)) {

    path <- file.path(
      outdir,
      paste0(
        label,
        "_db",
        i,
        "_unique.txt"
      )
    )

    writeLines(
      unique(mat_acc[, i]),
      con = path
    )

    paths[i] <- path
  }

  paths
}


.write_accession_counts <- function(
    db,
    indices,
    outdir,
    label
) {

  n_boot <- ncol(indices)
  n <- nrow(db)

  paths <- character(n_boot)

  for (i in seq_len(n_boot)) {

    counts <- tabulate(
      indices[, i],
      nbins = n
    )

    ids <- which(counts > 0L)

    out <- tibble::tibble(
      accession = db$accession[ids],
      name = db$name[ids],
      taxid = db$taxid[ids],
      filename = db$filename[ids],
      count = counts[ids]
    ) |>
      dplyr::arrange(
        dplyr::desc(count),
        accession
      )

    path <- file.path(
      outdir,
      paste0(
        label,
        "_db",
        i,
        "_accession_counts.tsv"
      )
    )

    readr::write_tsv(
      out,
      path
    )

    paths[i] <- path
  }

  paths
}


.write_taxid_counts <- function(
    mat_tax,
    outdir,
    label
) {

  n_boot <- ncol(mat_tax)
  paths <- character(n_boot)

  for (i in seq_len(n_boot)) {

    out <- tibble::tibble(
      taxid = mat_tax[, i]
    ) |>
      dplyr::count(
        taxid,
        name = "count"
      ) |>
      dplyr::arrange(
        dplyr::desc(count),
        taxid
      )

    path <- file.path(
      outdir,
      paste0(
        label,
        "_db",
        i,
        "_taxid_counts.tsv"
      )
    )

    readr::write_tsv(
      out,
      path
    )

    paths[i] <- path
  }

  paths
}


.write_all_bootstrap_accessions <- function(
    mat_acc,
    outdir,
    label
) {

  accessions <- sort(
    unique(c(mat_acc))
  )

  path <- file.path(
    outdir,
    paste0(
      label,
      "_all_accessions.txt"
    )
  )

  writeLines(
    accessions,
    con = path
  )

  list(
    path = path,
    accessions = accessions
  )
}


# Public function ---------------------------------------------------------

#' Generate bootstrap samples of genome assemblies
#'
#' Generates bootstrap samples from an assembly table produced by
#' [dedup_db()] or [filter_db()]. Each bootstrap replicate samples the
#' input assemblies with replacement and contains the same number of draws
#' as the original input table.
#'
#' @param db A data frame produced by [dedup_db()] or [filter_db()].
#' @param n_boot Number of bootstrap replicates to generate.
#' @param label Character label used to name bootstrap replicates and
#'   output files.
#' @param seed Optional integer random seed for reproducible sampling.
#' @param outdir Optional output directory. If `NULL`, defaults to
#'   `"<label>_boot"`.
#' @param outputs Character vector specifying output types. Supported
#'   values are `"composite"`, `"files"`, `"unique_files"`,
#'   `"accession_counts"`, `"taxid_counts"`, and `"all_accessions"`.
#' @param verbose Logical; print a summary of generated outputs.
#'
#' @return Invisibly returns a list containing run information and paths
#'   to generated output files.
#'
#' @details
#' `"files"` writes the raw bootstrap draws, including repeated accessions.
#' `"unique_files"` writes each accession only once per bootstrap replicate
#' and is suitable for identifying the genomes needed to construct each
#' database. `"accession_counts"` preserves the bootstrap multiplicity of
#' each sampled assembly.
#'
#' @examples
#' \dontrun{
#' db <- dedup_db("ncbi_genomes.tsv")
#'
#' db <- filter_db(
#'   db,
#'   assembly_level = "Complete Genome"
#' )
#'
#' boot <- bootstrap_db(
#'   db,
#'   n_boot = 100,
#'   label = "rickettsia",
#'   seed = 42
#' )
#' }
#'
#' @export
bootstrap_db <- function(
    db,
    n_boot,
    label,
    seed = NULL,
    outdir = NULL,
    outputs = c("composite", "files"),
    verbose = TRUE
) {

  db <- .validate_bootstrap_db(db)

  n_boot <- .validate_n_boot(n_boot)

  seed <- .validate_bootstrap_seed(seed)

  label_safe <- .sanitize_bootstrap_label(label)

  outputs <- .validate_bootstrap_outputs(
    outputs
  )

  if (is.null(outdir)) {
    outdir <- paste0(
      label_safe,
      "_boot"
    )
  }

  if (
    !is.character(outdir) ||
    length(outdir) != 1L ||
    is.na(outdir) ||
    stringr::str_trim(outdir) == ""
  ) {
    stop(
      "outdir must be a single non-empty path.",
      call. = FALSE
    )
  }

  outdir <- .make_bootstrap_outdir(
    outdir
  )

  boot <- .make_bootstrap_matrices(
    db = db,
    n_boot = n_boot,
    seed = seed
  )

  mat_acc <- .set_bootstrap_names(
    boot$accessions,
    label_safe
  )

  mat_tax <- .set_bootstrap_names(
    boot$taxids,
    label_safe
  )

  results <- list(
    label = label,
    label_safe = label_safe,
    outdir = outdir,
    n_assemblies = nrow(db),
    n_boot = n_boot,
    seed = seed,
    outputs = outputs
  )


  # Composite matrices ---------------------------------------------------

  if ("composite" %in% outputs) {

    paths <- .write_bootstrap_composites(
      mat_acc,
      mat_tax,
      outdir,
      label_safe
    )

    results$composite_accessions <-
      paths$accessions

    results$composite_taxids <-
      paths$taxids
  }


  # Raw bootstrap accession lists ---------------------------------------

  if ("files" %in% outputs) {

    results$bootstrap_files <-
      .write_bootstrap_files(
        mat_acc,
        outdir,
        label_safe
      )
  }


  # Unique accession lists ----------------------------------------------

  if ("unique_files" %in% outputs) {

    results$unique_bootstrap_files <-
      .write_unique_bootstrap_files(
        mat_acc,
        outdir,
        label_safe
      )
  }


  # Accession counts -----------------------------------------------------

  if ("accession_counts" %in% outputs) {

    results$accession_count_tables <-
      .write_accession_counts(
        db,
        boot$indices,
        outdir,
        label_safe
      )
  }


  # TaxID counts ---------------------------------------------------------

  if ("taxid_counts" %in% outputs) {

    results$taxid_count_tables <-
      .write_taxid_counts(
        mat_tax,
        outdir,
        label_safe
      )
  }


  # Union of accessions sampled anywhere --------------------------------

  if ("all_accessions" %in% outputs) {

    all_acc <- .write_all_bootstrap_accessions(
      mat_acc,
      outdir,
      label_safe
    )

    results$all_accessions_file <-
      all_acc$path

    results$all_accessions <-
      all_acc$accessions
  }


  # Report ---------------------------------------------------------------

  if (verbose) {

    cat(
      "\nBootstrap summary:\n",
      "  Input assemblies:     ", nrow(db), "\n",
      "  Bootstrap replicates: ", n_boot, "\n",
      "  Label:                ", label_safe, "\n",
      "  Output directory:     ", outdir, "\n",
      "  Outputs:              ", paste(outputs, collapse = ", "), "\n\n",
      sep = ""
    )
  }

  invisible(results)
}
