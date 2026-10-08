# Helper -----------------------------------------------------------------

# Write a small NCBI-style metadata table to a temporary TSV file.
# Each test can modify the table as needed without relying on external files.
write_test_metadata <- function(db) {
  path <- tempfile(fileext = ".tsv")

  readr::write_tsv(
    db,
    path,
    na = ""
  )

  path
}


# Standard test data ------------------------------------------------------

metadata <- tibble::tibble(
  assembly_name = c(
    "Assembly_A",
    "Assembly_A",
    "Assembly_B",
    "Assembly_C"
  ),
  assembly_accession = c(
    "GCA_000001.1",
    "GCF_000001.1",
    "GCA_000002.1",
    "GCF_000003.1"
  ),
  assembly_paired_assembly_accession = c(
    "GCF_000001.1",
    "GCA_000001.1",
    NA,
    NA
  ),
  organism_name = c(
    "Organism A",
    "Organism A",
    "Organism B",
    "Organism C"
  ),
  organism_taxonomic_id = c(
    101,
    101,
    102,
    103
  )
)


# Tests ------------------------------------------------------------------

test_that("dedup_db prefers RefSeq over paired GenBank assemblies", {

  infile <- write_test_metadata(metadata)

  result <- dedup_db(
    infile,
    report = FALSE
  )

  expect_true(
    "GCF_000001.1" %in% result$accession
  )

  expect_false(
    "GCA_000001.1" %in% result$accession
  )
})


test_that("dedup_db retains unpaired assemblies", {

  infile <- write_test_metadata(metadata)

  result <- dedup_db(
    infile,
    report = FALSE
  )

  expect_equal(
    sort(result$accession),
    sort(c(
      "GCF_000001.1",
      "GCA_000002.1",
      "GCF_000003.1"
    ))
  )

  expect_equal(
    nrow(result),
    3
  )
})


test_that("dedup_db adds standardized KrakPak columns", {

  infile <- write_test_metadata(metadata)

  result <- dedup_db(
    infile,
    report = FALSE
  )

  expect_true(
    all(
      c(
        "accession",
        "name",
        "taxid",
        "filename"
      ) %in% names(result)
    )
  )

  expect_equal(
    result$accession,
    result$assembly_accession
  )

  expect_equal(
    result$name,
    result$organism_name
  )

  expect_equal(
    result$taxid,
    result$organism_taxonomic_id
  )
})


test_that("dedup_db creates expected genome filenames", {

  infile <- write_test_metadata(metadata)

  result <- dedup_db(
    infile,
    report = FALSE
  )

  row_b <- result[
    result$accession == "GCA_000002.1",
  ]

  expect_equal(
    row_b$filename,
    "GCA_000002.1_Assembly_B_genomic.fna"
  )
})


test_that("dedup_db preserves additional metadata columns", {

  metadata_extra <- metadata

  metadata_extra$assembly_level <- c(
    "Complete Genome",
    "Complete Genome",
    "Chromosome",
    "Scaffold"
  )

  metadata_extra$contig_n50 <- c(
    100000,
    100000,
    50000,
    20000
  )

  infile <- write_test_metadata(metadata_extra)

  result <- dedup_db(
    infile,
    report = FALSE
  )

  expect_true(
    "assembly_level" %in% names(result)
  )

  expect_true(
    "contig_n50" %in% names(result)
  )
})


test_that("dedup_db errors when the input file does not exist", {

  expect_error(
    dedup_db(
      "this_file_does_not_exist.tsv",
      report = FALSE
    ),
    "Input file does not exist"
  )
})


test_that("dedup_db errors when infile is invalid", {

  expect_error(
    dedup_db(
      NULL,
      report = FALSE
    ),
    "infile must be a path"
  )

  expect_error(
    dedup_db(
      "",
      report = FALSE
    ),
    "infile must be a path"
  )
})


test_that("dedup_db errors when required columns are missing", {

  bad_metadata <- metadata |>
    dplyr::select(-organism_taxonomic_id)

  infile <- write_test_metadata(bad_metadata)

  expect_error(
    dedup_db(
      infile,
      report = FALSE
    ),
    "Missing required column"
  )
})


test_that("dedup_db errors when required values are missing", {

  bad_metadata <- metadata

  bad_metadata$organism_name[3] <- NA

  infile <- write_test_metadata(bad_metadata)

  expect_error(
    dedup_db(
      infile,
      report = FALSE
    ),
    "contains missing or empty values"
  )
})


test_that("dedup_db accepts missing paired assembly accessions", {

  unpaired <- tibble::tibble(
    assembly_name = c(
      "Assembly_A",
      "Assembly_B"
    ),
    assembly_accession = c(
      "GCA_000010.1",
      "GCF_000020.1"
    ),
    assembly_paired_assembly_accession = c(
      NA,
      NA
    ),
    organism_name = c(
      "Organism A",
      "Organism B"
    ),
    organism_taxonomic_id = c(
      201,
      202
    )
  )

  infile <- write_test_metadata(unpaired)

  result <- dedup_db(
    infile,
    report = FALSE
  )

  expect_equal(
    result$accession,
    c(
      "GCA_000010.1",
      "GCF_000020.1"
    )
  )
})
