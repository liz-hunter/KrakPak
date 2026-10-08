# KrakPak

**KrakPak** is an R package for building and evaluating bootstrap [Kraken2](https://github.com/DerrickWood/kraken2) databases.

KrakPak provides tools to prepare NCBI genome assembly metadata, generate reproducible bootstrap samples for database construction, read Kraken2 `inspect` reports, quantify variation in taxonomic representation across replicate databases, select taxa for in silico evaluation, and build classification manifests for downstream workflows.

> **Status:** KrakPak is a work in progress. 

## What KrakPak does

KrakPak currently supports six parts of a bootstrap Kraken database workflow:

1. **Prepare assembly metadata**
   - remove redundant paired GenBank/RefSeq assemblies with `dedup_db()`
   - filter assemblies using assembly level and quality metrics with `filter_db()`

2. **Generate bootstrap database inputs**
   - sample assemblies with replacement using `bootstrap_db()`
   - produce accession, taxid, and file manifests for downstream database construction

3. **Summarize bootstrap composition**
   - calculate taxid-based diversity statistics with `bootstrap_stats()`
   - visualize diversity across bootstrap replicates with `plot_bootstrap_stats()`

4. **Analyze Kraken2 inspect reports**
   - read and filter replicate reports with `read_inspect()`
   - extract descendants of a taxon with `inspect_children()`
   - reshape replicate data with `inspect_wide()`
   - quantify variability in taxon representation with `inspect_stats()`

5. **Select taxa for in silico evaluation**
   - classify taxa into low, intermediate, high, or other variability classes with `classify_inspect_taxa()`
   - select taxa by ranked or reproducible random sampling with `select_inspect_taxa()`
   - run both steps together with `make_insilico_taxa()`

6. **Build classification manifests**
   - match selected accessions to paired in silico FASTQ files with `make_manifest()`
   - expand each selected genome across replicate Kraken2 databases for downstream workflow execution

## Installation

KrakPak is currently available from GitHub.

```r
# install.packages("remotes")
remotes::install_github("liz-hunter/KrakPak")
```

For package development, clone the repository, open `KrakPak.Rproj` in RStudio, and load the development version with:

```r
devtools::load_all()
```

## Basic workflow

### 1. Prepare genome assembly metadata

`dedup_db()` reads an NCBI genome metadata TSV and removes redundant paired GenBank/RefSeq assemblies, retaining the RefSeq (`GCF_`) assembly when both versions are present.

```r
library(KrakPak)

db <- dedup_db("ncbi_genomes.tsv")
```

Optional assembly-quality filters can then be applied:

```r
db_filtered <- filter_db(
  db,
  assembly_level = "Complete Genome",
  min_checkm_completeness = 95,
  max_checkm_contamination = 5
)
```

Missing quality values are retained by default; use `keep_missing = FALSE` if assemblies lacking requested quality metrics should be removed.

### 2. Generate bootstrap replicates

`bootstrap_db()` samples assemblies with replacement. Each replicate contains the same number of draws as the input assembly table.

```r
boot <- bootstrap_db(
  db_filtered,
  n_boot = 100,
  label = "rickettsia",
  seed = 42
)
```

Several output formats are available, including composite bootstrap tables, per-replicate file lists, unique accession lists, accession counts, and taxid counts.

For example:

```r
boot <- bootstrap_db(
  db_filtered,
  n_boot = 100,
  label = "rickettsia",
  seed = 42,
  outputs = c(
    "composite",
    "unique_files",
    "accession_counts",
    "taxid_counts"
  )
)
```

### 3. Examine bootstrap diversity

When composite bootstrap output is generated, `bootstrap_stats()` calculates taxid-based diversity statistics for each replicate:

```r
stats <- bootstrap_stats(boot)

stats$summary
```

Reported metrics include:

- number of sampled assemblies
- number of unique taxids
- Shannon diversity
- Simpson dominance
- Gini-Simpson diversity
- inverse Simpson diversity
- Pielou evenness

The distributions can be visualized directly:

```r
plot_bootstrap_stats(stats)
```

or by replicate:

```r
plot_bootstrap_stats(
  stats,
  type = "replicate"
)
```

## Working with Kraken2 inspect reports

`read_inspect()` reads one or more Kraken2 inspect reports from a directory and combines them into a single tibble. Replicate/database identifiers are derived from the filenames.

```r
inspect <- read_inspect(
  directory = "path/to/inspect_reports",
  prefix = "rickettsia",
  level = "species"
)
```

KrakPak retains Kraken2's reported composition percentage and also calculates an exact percentage from the inclusive minimizer counts.

Taxonomic ranks can be selected with `level`, and specific taxids can be included or excluded:

```r
inspect <- read_inspect(
  directory = "path/to/inspect_reports",
  prefix = "rickettsia",
  level = "genus",
  exclude_taxids = c(12345, 67890)
)
```

### Inspect taxonomic descendants

To work with a complete Kraken taxonomy hierarchy, read all ranks including intermediate Kraken ranks:

```r
inspect_all <- read_inspect(
  directory = "path/to/inspect_reports",
  prefix = "rickettsia",
  level = "all",
  fuzzy = TRUE
)
```

Then extract a taxon's descendants:

```r
children <- inspect_children(
  inspect_all,
  taxid = 780
)
```

Use `immediate = TRUE` to return only immediate children.

### Compare replicate representation

`inspect_stats()` summarizes how consistently taxa are represented across replicate databases:

```r
inspect_summary <- inspect_stats(inspect)
```

The output includes replicate presence, mean and median representation, standard deviation, coefficient of variation, range, interquartile range, and inferential relative variance (InfRV).

To create a taxon-by-replicate table:

```r
inspect_matrix <- inspect_wide(
  inspect,
  value = "incl_min_count"
)
```

## Selecting taxa for in silico evaluation

KrakPak can use the output of `inspect_stats()` to define reproducible variability classes and select taxa for downstream in silico testing.

### Classify taxa by replicate variability

`classify_inspect_taxa()` assigns taxa to low, intermediate, high, or other variability classes using intersections of `report_rate`, coefficient of variation (`cv`), and an optional minimum mean representation.

The default criteria are:

- **Low:** `report_rate >= 0.9` and `cv <= 40`
- **Intermediate:** `report_rate` between `0.6` and `0.8` and `cv` between `50` and `100`
- **High:** `report_rate <= 0.3` and `cv >= 160`
- **Other:** taxa that do not satisfy one of the three designed classes

```r
classified <- classify_inspect_taxa(
  inspect_summary,
  tax_level = "S"
)
```

All thresholds are user-configurable. For example:

```r
classified <- classify_inspect_taxa(
  inspect_summary,
  tax_level = "S",
  low_report_min = 0.95,
  low_cv_max = 30,
  mid_report_min = 0.6,
  mid_report_max = 0.8,
  mid_cv_min = 50,
  mid_cv_max = 100,
  high_report_max = 0.2,
  high_cv_min = 200,
  min_mean = 0.01
)
```

### Select taxa within variability classes

`select_inspect_taxa()` selects a requested number of taxa from each class. Each class can be selected either by ranking or by reproducible random sampling.

By default:

- low variability taxa are selected by ranking, favoring better represented taxa
- intermediate variability taxa are selected randomly
- high variability taxa are selected by ranking, favoring better represented taxa within the unstable pool

```r
selected <- select_inspect_taxa(
  classified,
  n_per_group = 10,
  seed = 20260917,
  low_method = "ranked",
  intermediate_method = "random",
  high_method = "ranked"
)
```

Different sample sizes can be requested for each group:

```r
selected <- select_inspect_taxa(
  classified,
  n_per_group = c(
    low = 10,
    intermediate = 20,
    high = 10
  ),
  seed = 20260917
)
```

### Run classification and selection together

`make_insilico_taxa()` is a convenience wrapper that performs classification and selection in one step:

```r
selection <- make_insilico_taxa(
  inspect_summary,
  tax_level = "S",
  n_per_group = 10,
  seed = 20260917
)
```

The returned list contains:

- `classified`: all taxa at the requested taxonomic level with variability classes
- `candidates`: taxa belonging to the low, intermediate, or high candidate pools
- `selected`: taxa chosen for in silico evaluation
- `summary`: candidate and selected counts by variability class

```r
selection$selected
selection$summary
```

## Building a classification manifest

After in silico reads have been generated, `make_manifest()` matches selected genome accessions to paired FASTQ files in a designated directory and expands each genome across replicate Kraken2 databases.

FASTQ filenames must begin with an NCBI assembly accession and end in `_R1.fastq`, `_R2.fastq`, `_R1.fastq.gz`, or `_R2.fastq.gz`.

For example:

```text
GCF_000001.1_simulated_R1.fastq.gz
GCF_000001.1_simulated_R2.fastq.gz
```

Build a manifest directly from the taxa selected above:

```r
manifest <- make_manifest(
  selection$selected,
  reads_dir = "path/to/insilico_reads",
  n_db_reps = 10,
  out = "classification_manifest.tsv"
)
```

The resulting table contains the selected taxon metadata, matched `R1` and `R2` paths, and a `db_rep` column identifying the replicate Kraken2 database to use for each classification task.

This manifest is intended to provide a clean handoff from KrakPak's sample-selection logic to an external workflow manager such as Nextflow.

## Main functions

| Function | Purpose |
|---|---|
| `dedup_db()` | Remove redundant paired GenBank/RefSeq assemblies from NCBI metadata |
| `filter_db()` | Filter assemblies using assembly level and quality metrics |
| `bootstrap_db()` | Generate bootstrap samples and database-construction manifests |
| `bootstrap_stats()` | Calculate taxid-based diversity statistics across bootstrap replicates |
| `plot_bootstrap_stats()` | Plot bootstrap diversity statistics |
| `read_inspect()` | Read, combine, and filter Kraken2 inspect reports |
| `inspect_children()` | Extract descendants of a specified taxid from inspect reports |
| `inspect_wide()` | Reshape replicate inspect data to taxon-by-replicate wide format |
| `inspect_stats()` | Calculate variability statistics across inspect replicates |
| `classify_inspect_taxa()` | Classify taxa using report-rate and CV thresholds |
| `select_inspect_taxa()` | Select taxa from variability classes by ranked or random sampling |
| `make_insilico_taxa()` | Classify and select taxa for in silico evaluation in one step |
| `make_manifest()` | Match selected accessions to paired reads and expand classification tasks across DB replicates |

## Development

KrakPak is being developed as a research tool for evaluating how genome selection and database composition affect Kraken database behavior. Current development is focused on strengthening the R package interface, testing, documentation, reproducible sample selection, and downstream comparison utilities.

To check a local development copy:

```r
devtools::document()
devtools::test()
devtools::check()
```

Bug reports, suggestions, and contributions can be submitted through the GitHub repository.

## License

KrakPak is released under the MIT License.
