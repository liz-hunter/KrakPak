#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import re
import sys


# ============================================================
# Taxonomy loading
# ============================================================

def load_taxonomy(nodes_dmp):
    """
    Load parent and rank information from NCBI nodes.dmp.
    """

    parent = {}
    rank = {}

    with Path(nodes_dmp).open() as f:
        for line in f:
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
    Load old_taxid -> new_taxid mappings.
    """

    merged = {}

    path = Path(merged_dmp)

    if not path.exists():
        return merged

    with path.open() as f:
        for line in f:
            parts = line.split("|")

            if len(parts) < 2:
                continue

            old_taxid = parts[0].strip()
            new_taxid = parts[1].strip()

            merged[old_taxid] = new_taxid

    return merged


# ============================================================
# Taxid resolution
# ============================================================

def resolve_taxid(taxid, parent_map, merged_map):

    if pd.isna(taxid):
        return None

    taxid = str(taxid).strip()

    if taxid in ("", "0", "NA", "nan", "<NA>"):
        return None

    if taxid in parent_map:
        return taxid

    seen = set()

    while taxid in merged_map and taxid not in seen:

        seen.add(taxid)
        taxid = merged_map[taxid]

        if taxid in parent_map:
            return taxid

    return None


# ============================================================
# Project a taxid upward to requested rank
# ============================================================

def project_to_rank(
    taxid,
    target_rank,
    parent_map,
    rank_map,
    merged_map,
    cache,
):
    """
    Walk upward until target_rank is reached.

    If the taxid is already at target_rank, return it.

    If no ancestor at target_rank exists, return the resolved
    original taxid. This is important for calls ABOVE the target
    rank. Example:

        target = species
        prediction = genus

    The genus stays a genus, so its distance from the true species
    remains meaningful.
    """

    cache_key = (taxid, target_rank)

    if cache_key in cache:
        return cache[cache_key]

    resolved = resolve_taxid(
        taxid,
        parent_map,
        merged_map,
    )

    if resolved is None:
        cache[cache_key] = None
        return None

    cur = resolved
    seen = set()

    while cur in parent_map and cur not in seen:

        seen.add(cur)

        if rank_map.get(cur) == target_rank:
            cache[cache_key] = cur
            return cur

        par = parent_map[cur]

        if par == cur:
            break

        cur = par

    # No ancestor at requested rank.
    # Keep original resolved taxid.
    cache[cache_key] = resolved

    return resolved


# ============================================================
# Lineage + edge distance
# ============================================================

def lineage(taxid, parent_map, cache):

    if taxid in cache:
        return cache[taxid]

    lin = []
    cur = taxid
    seen = set()

    while cur in parent_map and cur not in seen:

        seen.add(cur)
        lin.append(cur)

        par = parent_map[cur]

        if par == cur:
            break

        cur = par

    cache[taxid] = lin

    return lin


def edge_distance(
    true_taxid,
    pred_taxid,
    target_rank,
    parent_map,
    rank_map,
    merged_map,
    rank_cache,
    lineage_cache,
):

    # --------------------------------------------------------
    # First normalize BOTH true and predicted taxids to the
    # requested analysis rank where possible.
    # --------------------------------------------------------

    true_norm = project_to_rank(
        true_taxid,
        target_rank,
        parent_map,
        rank_map,
        merged_map,
        rank_cache,
    )

    pred_norm = project_to_rank(
        pred_taxid,
        target_rank,
        parent_map,
        rank_map,
        merged_map,
        rank_cache,
    )

    if true_norm is None or pred_norm is None:
        return None, true_norm, pred_norm

    # Child classifications now collapse to same rank/node.
    if true_norm == pred_norm:
        return 0, true_norm, pred_norm

    true_lin = lineage(
        true_norm,
        parent_map,
        lineage_cache,
    )

    pred_lin = lineage(
        pred_norm,
        parent_map,
        lineage_cache,
    )

    pred_positions = {
        taxid: i
        for i, taxid in enumerate(pred_lin)
    }

    for i, ancestor in enumerate(true_lin):

        j = pred_positions.get(ancestor)

        if j is not None:
            return i + j, true_norm, pred_norm

    return None, true_norm, pred_norm


# ============================================================
# Main
# ============================================================

def main():

    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--input-tsv",
        required=True,
        help="Wide per-read Kraken summary TSV",
    )

    ap.add_argument(
        "--key",
        required=True,
        help="extremes_key_simple.tsv",
    )

    ap.add_argument(
        "--myfiles",
        required=True,
        help="Value of $MYFILES",
    )

    ap.add_argument(
        "--rank",
        default="species",
        help="Rank at which classifications are evaluated "
             "(default: species)",
    )

    ap.add_argument(
        "--output-dir",
        required=True,
    )

    args = ap.parse_args()

    input_file = Path(args.input_tsv)
    output_dir = Path(args.output_dir)
    myfiles = Path(args.myfiles)

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    per_read_dir = output_dir / "per_read"
    summary_dir = output_dir / "summary_parts"

    per_read_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    summary_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ========================================================
    # Read key
    # ========================================================

    key = pd.read_csv(
        args.key,
        sep="\t",
        dtype={
            "dataset": "string",
            "taxid": "string",
            "accession": "string",
        },
    )

    # ========================================================
    # Determine accession from filename
    #
    # Files created by your current summarizer are:
    # ACCESSION_dataset_classification.tsv
    # ========================================================

    matching = key[
        key["accession"].apply(
            lambda x: input_file.name.startswith(str(x))
        )
    ]

    if len(matching) != 1:
        sys.exit(
            f"ERROR: expected exactly one accession from the key "
            f"to match {input_file.name}; found {len(matching)}"
        )

    key_row = matching.iloc[0]

    accession = str(key_row["accession"])
    dataset = str(key_row["dataset"])
    true_taxid = str(key_row["taxid"])

    print("=" * 70)
    print(f"Input:       {input_file.name}")
    print(f"Accession:   {accession}")
    print(f"Dataset:     {dataset}")
    print(f"True taxid:  {true_taxid}")
    print(f"Target rank: {args.rank}")
    print("=" * 70)

    # ========================================================
    # Select full taxonomy
    # ========================================================

    taxonomy_dirs = {
        "entero":
            myfiles
            / "entropy_entero"
            / "entero_full"
            / "taxonomy",

        "stan":
            myfiles
            / "entropy_stan"
            / "stan_full"
            / "taxonomy",

        "virid":
            myfiles
            / "entropy_virid"
            / "virid_full"
            / "taxonomy",
    }

    if dataset not in taxonomy_dirs:
        sys.exit(
            f"ERROR: unknown dataset '{dataset}'"
        )

    taxonomy_dir = taxonomy_dirs[dataset]

    nodes_dmp = taxonomy_dir / "nodes.dmp"
    merged_dmp = taxonomy_dir / "merged.dmp"

    if not nodes_dmp.exists():
        sys.exit(
            f"ERROR: nodes.dmp not found: {nodes_dmp}"
        )

    # ========================================================
    # Load taxonomy
    # ========================================================

    print(f"Loading taxonomy: {nodes_dmp}")

    parent_map, rank_map = load_taxonomy(
        nodes_dmp
    )

    merged_map = load_merged(
        merged_dmp
    )

    print(
        f"Loaded {len(parent_map):,} nodes "
        f"and {len(merged_map):,} merged taxids"
    )

    rank_cache = {}
    lineage_cache = {}

    # ========================================================
    # Normalize TRUE taxid once
    # ========================================================

    normalized_true = project_to_rank(
        true_taxid,
        args.rank,
        parent_map,
        rank_map,
        merged_map,
        rank_cache,
    )

    if normalized_true is None:
        sys.exit(
            f"ERROR: true taxid {true_taxid} cannot be resolved"
        )

    print(
        f"True taxid at {args.rank}: "
        f"{true_taxid} -> {normalized_true}"
    )

    # ========================================================
    # Read wide Kraken summary
    # ========================================================

    df = pd.read_csv(
        input_file,
        sep="\t",
        dtype="string",
    )

    # ========================================================
    # Find db columns
    # ========================================================

    pattern = re.compile(
        rf"^{re.escape(dataset)}_db(\d+)$"
    )

    db_columns = []

    for col in df.columns:

        m = pattern.match(col)

        if m:
            db_columns.append(
                (int(m.group(1)), col)
            )

    db_columns.sort()

    if not db_columns:
        sys.exit(
            f"ERROR: no {dataset}_dbN columns found"
        )

    summary_rows = []

    # ========================================================
    # Process each Kraken database
    # ========================================================

    for db_number, pred_col in db_columns:

        print(f"\nProcessing {pred_col}")

        edge_col = (
            f"{pred_col}_edge_distance"
        )

        # ----------------------------------------------------
        # Compute only once per unique predicted taxid
        # ----------------------------------------------------

        unique_predictions = (
            df[pred_col]
            .dropna()
            .unique()
        )

        distance_lookup = {}
        normalized_lookup = {}

        for pred_taxid in unique_predictions:

            pred_taxid = str(pred_taxid)

            dist, _, pred_norm = edge_distance(
                true_taxid=true_taxid,
                pred_taxid=pred_taxid,
                target_rank=args.rank,
                parent_map=parent_map,
                rank_map=rank_map,
                merged_map=merged_map,
                rank_cache=rank_cache,
                lineage_cache=lineage_cache,
            )

            distance_lookup[pred_taxid] = dist
            normalized_lookup[pred_taxid] = pred_norm

        # ----------------------------------------------------
        # Map distances back to reads
        # ----------------------------------------------------

        df[edge_col] = (
            df[pred_col]
            .map(distance_lookup)
            .astype("Int64")
        )

        # ----------------------------------------------------
        # Summary
        # ----------------------------------------------------

        n_reads = len(df)

        unclassified = (
            df[pred_col].isna()
            | (df[pred_col] == "0")
        )

        n_unclassified = int(
            unclassified.sum()
        )

        n_valid = int(
            df[edge_col].notna().sum()
        )

        n_correct = int(
            (df[edge_col] == 0).sum()
        )

        n_no_distance = (
            n_reads
            - n_unclassified
            - n_valid
        )

        mean_distance = (
            float(df[edge_col].mean())
            if n_valid > 0
            else float("nan")
        )

        pct_correct = (
            100.0 * n_correct / n_reads
            if n_reads > 0
            else float("nan")
        )

        print(f"  reads:            {n_reads:,}")
        print(f"  valid distances:  {n_valid:,}")
        print(f"  correct @ rank:   {n_correct:,}")
        print(f"  unclassified:     {n_unclassified:,}")
        print(f"  no distance:      {n_no_distance:,}")
        print(f"  mean edges:       {mean_distance:.6f}")
        print(f"  percent correct:  {pct_correct:.4f}")

        summary_rows.append(
            {
                "taxid": true_taxid,
                "rank_taxid": normalized_true,
                "rank": args.rank,
                "db": db_number,
                "n_reads": n_reads,
                "n_valid": n_valid,
                "n_correct": n_correct,
                "n_unclassified": n_unclassified,
                "n_no_distance": n_no_distance,
                "pct_correct": pct_correct,
                "mean_edge_distance": mean_distance,
            }
        )

    # ========================================================
    # Write outputs
    # ========================================================

    stem = input_file.stem

    per_read_file = (
        per_read_dir
        / f"{stem}_{args.rank}_edges.tsv"
    )

    summary_file = (
        summary_dir
        / f"{accession}_{dataset}_{args.rank}_summary.tsv"
    )

    df.to_csv(
        per_read_file,
        sep="\t",
        index=False,
    )

    pd.DataFrame(
        summary_rows
    ).to_csv(
        summary_file,
        sep="\t",
        index=False,
    )

    print("\n" + "=" * 70)
    print(f"Wrote: {per_read_file}")
    print(f"Wrote: {summary_file}")
    print("=" * 70)


if __name__ == "__main__":
    main()
