#' Build an in silico classification manifest
#'
#' Matches selected taxa to paired in silico FASTQ files in a directory
#' and expands each sample across replicate Kraken databases.
#'
#' FASTQ filenames must begin with an NCBI assembly accession and end in
#' `_R1.fastq`, `_R2.fastq`, `_R1.fastq.gz`, or `_R2.fastq.gz`.
#'
#' @param selected A data frame of selected taxa. Must contain an
#'   `accession` column.
#' @param reads_dir Directory containing paired in silico FASTQ files.
#' @param n_db_reps Number of Kraken database replicates. Defaults to 10.
#' @param out Optional output TSV path. If `NULL`, no file is written.
#'
#' @return A tibble containing the selected taxa, matched `R1` and `R2`
#'   paths, and `db_rep`.
#'
#' @export
make_manifest <- function(
    selected,
    reads_dir,
    n_db_reps = 10,
    out = NULL
) {

  if (!is.data.frame(selected)) {
    stop(
      "selected must be a data frame.",
      call. = FALSE
    )
  }

  if (!"accession" %in% names(selected)) {
    stop(
      "selected must contain an 'accession' column.",
      call. = FALSE
    )
  }

  if (
    length(reads_dir) != 1L ||
    !is.character(reads_dir) ||
    is.na(reads_dir) ||
    !nzchar(reads_dir) ||
    !dir.exists(reads_dir)
  ) {
    stop(
      "reads_dir must be an existing directory.",
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

  fastq_files <- list.files(
    reads_dir,
    pattern = "\\.fastq(?:\\.gz)?$",
    full.names = TRUE
  )

  if (length(fastq_files) == 0L) {
    stop(
      "No FASTQ files found in reads_dir.",
      call. = FALSE
    )
  }

  fastq_names <- basename(fastq_files)

  pattern <- paste0(
    "^(GC[AF]_[0-9]+\\.[0-9]+)",
    ".*_R([12])\\.fastq(?:\\.gz)?$"
  )

  matched <- stringr::str_match(
    fastq_names,
    pattern
  )

  keep <- !is.na(matched[, 1])

  reads <- tibble::tibble(
    accession = matched[keep, 2],
    mate = matched[keep, 3],
    path = fastq_files[keep]
  )

  if (nrow(reads) == 0L) {
    stop(
      "No FASTQ filenames matched the expected accession/R1/R2 pattern.",
      call. = FALSE
    )
  }

  duplicates <- reads |>
    dplyr::count(
      .data$accession,
      .data$mate,
      name = "n"
    ) |>
    dplyr::filter(
      .data$n > 1L
    )

  if (nrow(duplicates) > 0L) {

    first <- duplicates[1, ]

    stop(
      "Duplicate R",
      first$mate,
      " file for ",
      first$accession,
      ".",
      call. = FALSE
    )
  }

  reads <- reads |>
    tidyr::pivot_wider(
      names_from = "mate",
      values_from = "path",
      names_prefix = "R"
    )

  selected <- tibble::as_tibble(selected)

  manifest <- selected |>
    dplyr::left_join(
      reads,
      by = "accession"
    )

  missing_both <- manifest$accession[
    is.na(manifest$R1) &
      is.na(manifest$R2)
  ]

  if (length(missing_both) > 0L) {
    stop(
      "No FASTQ files found for accession(s): ",
      paste(
        missing_both,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  missing_r1 <- manifest$accession[
    is.na(manifest$R1)
  ]

  if (length(missing_r1) > 0L) {
    stop(
      "No R1 found for accession(s): ",
      paste(
        missing_r1,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  missing_r2 <- manifest$accession[
    is.na(manifest$R2)
  ]

  if (length(missing_r2) > 0L) {
    stop(
      "No R2 found for accession(s): ",
      paste(
        missing_r2,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  manifest <- manifest[
    rep(
      seq_len(nrow(manifest)),
      each = n_db_reps
    ),
    ,
    drop = FALSE
  ]

  manifest$db_rep <- rep(
    seq_len(n_db_reps),
    times = nrow(selected)
  )

  manifest <- tibble::as_tibble(manifest)

  if (!is.null(out)) {

    readr::write_tsv(
      manifest,
      out
    )
  }

  manifest
}
