from pathlib import Path
import re

import pandas as pd


# -------------------------------------------------------------------------
# Taxonomy loading
# -------------------------------------------------------------------------

def load_taxonomy(nodes_dmp):
    """
    Load parent and rank information from an NCBI nodes.dmp file.

    Parameters
    ----------
    nodes_dmp : str or pathlib.Path
        Path to nodes.dmp.

    Returns
    -------
    tuple
        Two dictionaries:

        parent
            taxid -> parent taxid

        rank
            taxid -> taxonomic rank
    """

    nodes_dmp = Path(nodes_dmp)

    if not nodes_dmp.exists():
        raise FileNotFoundError(
            f"nodes.dmp not found: {nodes_dmp}"
        )

    parent = {}
    rank = {}

    with nodes_dmp.open() as fh:
        for line in fh:

            parts = line.split("|")

            if len(parts) < 3:
                continue

            taxid = parts[0].strip()
            parent_taxid = parts[1].strip()
            tax_rank = parts[2].strip()

            parent[taxid] = parent_taxid
            rank[taxid] = tax_rank

    return parent, rank


def load_merged(merged_dmp):
    """
    Load old-taxid to current-taxid mappings from merged.dmp.

    If merged.dmp does not exist, an empty mapping is returned.

    Parameters
    ----------
    merged_dmp : str or pathlib.Path
        Path to merged.dmp.

    Returns
    -------
    dict
        old taxid -> current taxid
    """

    merged_dmp = Path(merged_dmp)

    merged = {}

    if not merged_dmp.exists():
        return merged

    with merged_dmp.open() as fh:
        for line in fh:

            parts = line.split("|")

            if len(parts) < 2:
                continue

            old_taxid = parts[0].strip()
            new_taxid = parts[1].strip()

            merged[old_taxid] = new_taxid

    return merged


# -------------------------------------------------------------------------
# Taxid handling
# -------------------------------------------------------------------------

def resolve_taxid(
    taxid,
    parent_map,
    merged_map,
):
    """
    Resolve a taxid against the current NCBI taxonomy.

    Current taxids are returned unchanged. Obsolete taxids are followed
    through merged.dmp until a current taxid is found.

    Unclassified, missing, zero, or unresolvable taxids return None.
    """

    if pd.isna(taxid):
        return None

    taxid = str(taxid).strip()

    if taxid in (
        "",
        "0",
        "NA",
        "nan",
        "<NA>",
        "unclassified",
    ):
        return None

    if taxid in parent_map:
        return taxid

    seen = set()

    while (
        taxid in merged_map
        and taxid not in seen
    ):

        seen.add(taxid)

        taxid = merged_map[taxid]

        if taxid in parent_map:
            return taxid

    return None


def project_to_rank(
    taxid,
    target_rank,
    parent_map,
    rank_map,
    merged_map,
    cache=None,
):
    """
    Project a taxid upward to a requested taxonomic rank.

    If the resolved taxid is already at the target rank, it is returned.

    If the taxid is below the target rank, its ancestor at the target rank
    is returned.

    If no ancestor at the target rank exists, the resolved original taxid
    is retained. This preserves meaningful distances for classifications
    made above the target rank.

    For example, when evaluating at species rank, a prediction made at
    genus rank remains a genus rather than becoming missing.
    """

    if cache is None:
        cache = {}

    cache_key = (
        str(taxid),
        target_rank,
    )

    if cache_key in cache:
        return cache[cache_key]

    resolved = resolve_taxid(
        taxid=taxid,
        parent_map=parent_map,
        merged_map=merged_map,
    )

    if resolved is None:
        cache[cache_key] = None
        return None

    current = resolved
    seen = set()

    while (
        current in parent_map
        and current not in seen
    ):

        seen.add(current)

        if rank_map.get(current) == target_rank:
            cache[cache_key] = current
            return current

        parent = parent_map[current]

        if parent == current:
            break

        current = parent

    # No ancestor at requested rank.
    # Keep the resolved original taxid.
    cache[cache_key] = resolved

    return resolved


# -------------------------------------------------------------------------
# Lineage and edge distance
# -------------------------------------------------------------------------

def lineage(
    taxid,
    parent_map,
    cache=None,
):
    """
    Return lineage from taxid upward toward the taxonomy root.

    The starting taxid is included as the first element.
    """

    if cache is None:
        cache = {}

    if taxid in cache:
        return cache[taxid]

    result = []

    current = taxid
    seen = set()

    while (
        current in parent_map
        and current not in seen
    ):

        seen.add(current)

        result.append(current)

        parent = parent_map[current]

        if parent == current:
            break

        current = parent

    cache[taxid] = result

    return result


def edge_distance(
    true_taxid,
    pred_taxid,
    target_rank,
    parent_map,
    rank_map,
    merged_map,
    rank_cache=None,
    lineage_cache=None,
):
    """
    Calculate taxonomy-tree edge distance between true and predicted taxids.

    Both taxids are projected to the requested analysis rank where possible.

    Parameters
    ----------
    true_taxid : str
        True taxonomic ID.

    pred_taxid : str
        Predicted taxonomic ID.

    target_rank : str
        Taxonomic rank at which the comparison should be evaluated.

    parent_map : dict
        taxid -> parent taxid

    rank_map : dict
        taxid -> rank

    merged_map : dict
        old taxid -> current taxid

    rank_cache : dict, optional
        Cache for rank projection.

    lineage_cache : dict, optional
        Cache for lineage calculations.

    Returns
    -------
    tuple
        distance, normalized true taxid, normalized predicted taxid

        distance is None when either taxid cannot be resolved or when the
        two lineages have no common node in the supplied taxonomy.
    """

    if rank_cache is None:
        rank_cache = {}

    if lineage_cache is None:
        lineage_cache = {}

    true_norm = project_to_rank(
        taxid=true_taxid,
        target_rank=target_rank,
        parent_map=parent_map,
        rank_map=rank_map,
        merged_map=merged_map,
        cache=rank_cache,
    )

    pred_norm = project_to_rank(
        taxid=pred_taxid,
        target_rank=target_rank,
        parent_map=parent_map,
        rank_map=rank_map,
        merged_map=merged_map,
        cache=rank_cache,
    )

    if (
        true_norm is None
        or pred_norm is None
    ):
        return (
            None,
            true_norm,
            pred_norm,
        )

    # Same taxon after projection to target rank.
    if true_norm == pred_norm:
        return (
            0,
            true_norm,
            pred_norm,
        )

    true_lineage = lineage(
        taxid=true_norm,
        parent_map=parent_map,
        cache=lineage_cache,
    )

    pred_lineage = lineage(
        taxid=pred_norm,
        parent_map=parent_map,
        cache=lineage_cache,
    )

    pred_positions = {
        taxid: index
        for index, taxid in enumerate(
            pred_lineage
        )
    }

    for true_index, ancestor in enumerate(
        true_lineage
    ):

        pred_index = pred_positions.get(
            ancestor
        )

        if pred_index is not None:

            distance = (
                true_index
                + pred_index
            )

            return (
                distance,
                true_norm,
                pred_norm,
            )

    return (
        None,
        true_norm,
        pred_norm,
    )


# -------------------------------------------------------------------------
# Classification-table helpers
# -------------------------------------------------------------------------

def _read_classification(
    classification,
):
    """
    Read or copy a KrakPak per-read classification table.
    """

    if isinstance(
        classification,
        pd.DataFrame,
    ):
        return classification.copy()

    path = Path(classification)

    if not path.exists():
        raise FileNotFoundError(
            f"Classification file not found: {path}"
        )

    return pd.read_csv(
        path,
        sep="\t",
        dtype="string",
    )


def _find_db_columns(
    df,
    dataset=None,
):
    """
    Identify columns named <dataset>_dbN.

    Returns
    -------
    tuple
        inferred/validated dataset name,
        list of (database number, column name)
    """

    if dataset is None:

        pattern = re.compile(
            r"^(.+)_db(\d+)$"
        )

        matches = []

        for column in df.columns:

            match = pattern.match(
                column
            )

            if match:
                matches.append(
                    (
                        match.group(1),
                        int(match.group(2)),
                        column,
                    )
                )

        if not matches:
            raise ValueError(
                "No classification columns matching "
                "'<dataset>_dbN' were found."
            )

        datasets = {
            item[0]
            for item in matches
        }

        if len(datasets) != 1:
            raise ValueError(
                "More than one dataset prefix was found in "
                "classification columns: "
                + ", ".join(
                    sorted(datasets)
                )
            )

        dataset = next(
            iter(datasets)
        )

        db_columns = [
            (
                db_number,
                column,
            )
            for prefix, db_number, column
            in matches
            if prefix == dataset
        ]

    else:

        dataset = str(
            dataset
        )

        pattern = re.compile(
            rf"^{re.escape(dataset)}_db(\d+)$"
        )

        db_columns = []

        for column in df.columns:

            match = pattern.match(
                column
            )

            if match:
                db_columns.append(
                    (
                        int(
                            match.group(1)
                        ),
                        column,
                    )
                )

    db_columns.sort(
        key=lambda x: x[0]
    )

    if not db_columns:
        raise ValueError(
            f"No {dataset}_dbN columns found."
        )

    return (
        dataset,
        db_columns,
    )


# -------------------------------------------------------------------------
# Main reusable function
# -------------------------------------------------------------------------

def edge_walk_classification(
    classification,
    taxonomy_dir,
    target_rank="species",
    dataset=None,
    output_dir=None,
):
    """
    Calculate taxonomic edge distances for a KrakPak classification table.

    Parameters
    ----------
    classification : str, pathlib.Path, or pandas.DataFrame
        Per-read classification table produced by KrakPak's
        summarize_classification().

        Expected columns include:

            read_id
            true_taxid
            <dataset>_db1
            <dataset>_db2
            ...

        Unclassified reads should contain the literal value
        "unclassified".

    taxonomy_dir : str or pathlib.Path
        Directory containing NCBI nodes.dmp and optionally merged.dmp.

    target_rank : str, default "species"
        Taxonomic rank at which classifications should be evaluated.

    dataset : str, optional
        Dataset prefix used in classification columns.

        If None, the prefix is inferred from columns named
        <dataset>_dbN.

    output_dir : str or pathlib.Path, optional
        If supplied, write per-read and summary TSV files here.

    Returns
    -------
    dict
        {
            "per_read": pandas.DataFrame,
            "summary": pandas.DataFrame,
            "dataset": str,
            "target_rank": str,
            "normalized_true_taxid": str,
            "per_read_file": pathlib.Path or None,
            "summary_file": pathlib.Path or None,
        }
    """

    taxonomy_dir = Path(
        taxonomy_dir
    )

    nodes_dmp = (
        taxonomy_dir
        / "nodes.dmp"
    )

    merged_dmp = (
        taxonomy_dir
        / "merged.dmp"
    )

    if not nodes_dmp.exists():
        raise FileNotFoundError(
            f"nodes.dmp not found: {nodes_dmp}"
        )

    # ------------------------------------------------------------------
    # Load taxonomy
    # ------------------------------------------------------------------

    parent_map, rank_map = (
        load_taxonomy(
            nodes_dmp
        )
    )

    merged_map = load_merged(
        merged_dmp
    )

    rank_cache = {}
    lineage_cache = {}

    # ------------------------------------------------------------------
    # Read classification table
    # ------------------------------------------------------------------

    df = _read_classification(
        classification
    )

    required = {
        "read_id",
        "true_taxid",
    }

    missing = (
        required
        - set(df.columns)
    )

    if missing:
        raise ValueError(
            "Classification table is missing required "
            "column(s): "
            + ", ".join(
                sorted(missing)
            )
        )

    # Require one true taxid per input table.
    true_taxids = (
        df["true_taxid"]
        .dropna()
        .astype(str)
        .unique()
    )

    if len(true_taxids) != 1:
        raise ValueError(
            "Classification table must contain exactly "
            "one unique true_taxid."
        )

    true_taxid = str(
        true_taxids[0]
    )

    # ------------------------------------------------------------------
    # Identify database columns
    # ------------------------------------------------------------------

    dataset, db_columns = (
        _find_db_columns(
            df=df,
            dataset=dataset,
        )
    )

    # ------------------------------------------------------------------
    # Normalize true taxid once
    # ------------------------------------------------------------------

    normalized_true = project_to_rank(
        taxid=true_taxid,
        target_rank=target_rank,
        parent_map=parent_map,
        rank_map=rank_map,
        merged_map=merged_map,
        cache=rank_cache,
    )

    if normalized_true is None:
        raise ValueError(
            f"True taxid {true_taxid} could not be "
            "resolved in the supplied taxonomy."
        )

    summary_rows = []

    # ------------------------------------------------------------------
    # Process database replicates
    # ------------------------------------------------------------------

    for db_number, pred_col in db_columns:

        edge_col = (
            f"{pred_col}_edge_distance"
        )

        predictions = (
            df[pred_col]
            .astype("string")
        )

        # Explicit Kraken outcome.
        unclassified = (
            predictions
            == "unclassified"
        )

        # Missing means the read is absent from this replicate's
        # classification table, which is different from explicitly
        # unclassified.
        missing_prediction = (
            predictions.isna()
        )

        # --------------------------------------------------------------
        # Compute each unique predicted taxid only once
        # --------------------------------------------------------------

        unique_predictions = (
            predictions[
                ~unclassified
                & ~missing_prediction
            ]
            .dropna()
            .unique()
        )

        distance_lookup = {}
        normalized_lookup = {}

        for pred_taxid in unique_predictions:

            pred_taxid = str(
                pred_taxid
            )

            (
                distance,
                _,
                pred_norm,
            ) = edge_distance(
                true_taxid=true_taxid,
                pred_taxid=pred_taxid,
                target_rank=target_rank,
                parent_map=parent_map,
                rank_map=rank_map,
                merged_map=merged_map,
                rank_cache=rank_cache,
                lineage_cache=lineage_cache,
            )

            distance_lookup[
                pred_taxid
            ] = distance

            normalized_lookup[
                pred_taxid
            ] = pred_norm

        # --------------------------------------------------------------
        # Map back to reads
        # --------------------------------------------------------------

        edge_values = predictions.map(
            distance_lookup
        )

        df[edge_col] = pd.array(
            edge_values,
            dtype="Int64",
        )

        # --------------------------------------------------------------
        # Summary
        # --------------------------------------------------------------

        n_reads = len(df)

        n_unclassified = int(
            unclassified.sum()
        )

        n_missing_prediction = int(
            missing_prediction.sum()
        )

        n_valid_distance = int(
            df[edge_col]
            .notna()
            .sum()
        )

        n_correct = int(
            (
                df[edge_col]
                == 0
            )
            .fillna(False)
            .sum()
        )

        n_incorrect_distance = int(
            (
                df[edge_col]
                > 0
            )
            .fillna(False)
            .sum()
        )

        # Classified taxids that could not be resolved in the
        # supplied taxonomy.
        classified_taxid = (
            ~unclassified
            & ~missing_prediction
        )

        no_distance = (
            classified_taxid
            & df[edge_col].isna()
        )

        n_no_distance = int(
            no_distance.sum()
        )

        mean_edge_distance = (
            float(
                df[edge_col].mean()
            )
            if n_valid_distance > 0
            else float("nan")
        )

        pct_correct_all = (
            100.0
            * n_correct
            / n_reads
            if n_reads > 0
            else float("nan")
        )

        pct_correct_valid = (
            100.0
            * n_correct
            / n_valid_distance
            if n_valid_distance > 0
            else float("nan")
        )

        summary_rows.append(
            {
                "true_taxid":
                    true_taxid,

                "rank_taxid":
                    normalized_true,

                "rank":
                    target_rank,

                "dataset":
                    dataset,

                "db":
                    db_number,

                "n_reads":
                    n_reads,

                "n_valid_distance":
                    n_valid_distance,

                "n_correct":
                    n_correct,

                "n_incorrect_distance":
                    n_incorrect_distance,

                "n_unclassified":
                    n_unclassified,

                "n_no_distance":
                    n_no_distance,

                "n_missing_prediction":
                    n_missing_prediction,

                "pct_correct_all":
                    pct_correct_all,

                "pct_correct_valid":
                    pct_correct_valid,

                "mean_edge_distance":
                    mean_edge_distance,
            }
        )

    summary = pd.DataFrame(
        summary_rows
    )

    # ------------------------------------------------------------------
    # Optional output
    # ------------------------------------------------------------------

    per_read_file = None
    summary_file = None

    if output_dir is not None:

        output_dir = Path(
            output_dir
        )

        output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        per_read_file = (
            output_dir
            / (
                f"{dataset}_"
                f"{true_taxid}_"
                f"{target_rank}_edges.tsv"
            )
        )

        summary_file = (
            output_dir
            / (
                f"{dataset}_"
                f"{true_taxid}_"
                f"{target_rank}_edge_summary.tsv"
            )
        )

        df.to_csv(
            per_read_file,
            sep="\t",
            index=False,
        )

        summary.to_csv(
            summary_file,
            sep="\t",
            index=False,
        )

    return {
        "per_read":
            df,

        "summary":
            summary,

        "dataset":
            dataset,

        "target_rank":
            target_rank,

        "normalized_true_taxid":
            normalized_true,

        "per_read_file":
            per_read_file,

        "summary_file":
            summary_file,
    }
