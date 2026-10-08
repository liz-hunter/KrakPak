#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import re
import sys


def extract_taxid(value):
    """
    Convert Kraken's taxon field to numeric taxid only.

    Handles both:
        4146
    and:
        Olea europaea (taxid 4146)
    """

    if pd.isna(value):
        return pd.NA

    value = str(value).strip()

    # Already just a numeric taxid
    if value.isdigit():
        return value

    # Kraken --use-names style:
    match = re.search(r"\(taxid\s+(\d+)\)\s*$", value)

    if match:
        return match.group(1)

    return pd.NA


def read_kraken_file(path):
    """
    Read only the Kraken columns we actually need:
      column 0 = classification status
      column 1 = read ID
      column 2 = predicted taxon

    The predicted taxon is normalized immediately to numeric taxid only.
    """

    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["status", "read_id", "pred_taxid"],
        dtype={
            "status": "string",
            "read_id": "string",
            "pred_taxid": "string",
        },
    )

    # Convert taxon names such as Olea europaea (taxid 4146)" to "4146"
    df["pred_taxid"] = (
        df["pred_taxid"]
        .map(extract_taxid)
        .astype("string")
    )

    # A classified read should always have a recognizable taxid
    bad = (
        (df["status"] == "C")
        & df["pred_taxid"].isna()
    )

    if bad.any():
        raise ValueError(
            f"{int(bad.sum()):,} classified reads in {path} "
            f"did not contain a recognizable taxid"
        )

    return df


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--key",
        required=True,
        help="Simplified TSV key with dataset, taxid, accession",
    )

    parser.add_argument(
        "--kraken-dir",
        required=True,
        help="Directory containing Kraken output files",
    )

    parser.add_argument(
        "--out-dir",
        required=True,
        help="Output directory",
    )

    parser.add_argument(
        "--task-id",
        type=int,
        required=True,
        help="Zero-based row number from key",
    )

    parser.add_argument(
        "--n-db",
        type=int,
        default=10,
        help="Number of Kraken database replicates (default: 10)",
    )

    args = parser.parse_args()

    key_file = Path(args.key)
    kraken_dir = Path(args.kraken_dir)
    out_dir = Path(args.out_dir)

    out_dir.mkdir(parents=True, exist_ok=True)

    summary_dir = out_dir / "summary_parts"
    summary_dir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------
    # Read key
    # --------------------------------------------------------

    key = pd.read_csv(
        key_file,
        sep="\t",
        dtype={
            "dataset": "string",
            "taxid": "string",
            "accession": "string",
        },
    )

    if args.task_id < 0 or args.task_id >= len(key):
        sys.exit(
            f"ERROR: task ID {args.task_id} is outside key range "
            f"0-{len(key) - 1}"
        )

    row = key.iloc[args.task_id]

    dataset = row["dataset"]
    true_taxid = row["taxid"]
    accession = row["accession"]

    print("=" * 70)
    print(f"Task ID:      {args.task_id}")
    print(f"Accession:    {accession}")
    print(f"Dataset:      {dataset}")
    print(f"True taxid:   {true_taxid}")
    print("=" * 70)

    # --------------------------------------------------------
    # Process db1-dbN
    # --------------------------------------------------------

    result = None
    summary_rows = []

    for db in range(1, args.n_db + 1):

        filename = f"{accession}_{dataset}_db{db}.txt"
        path = kraken_dir / filename

        if not path.exists():
            print(f"WARNING: missing file: {path}")

            summary_rows.append(
                {
                    "accession": accession,
                    "dataset": dataset,
                    "db": db,
                    "true_taxid": true_taxid,
                    "n_reads": pd.NA,
                    "n_classified": pd.NA,
                    "n_correct": pd.NA,
                    "n_incorrect": pd.NA,
                    "pct_correct": pd.NA,
                    "file_found": False,
                }
            )

            continue

        print(f"Reading db{db}: {filename}")

        x = read_kraken_file(path)

        # ----------------------------------------------------
        # Classification QC
        # ----------------------------------------------------

        n_reads = len(x)

        n_classified = int(
            (x["status"] == "C").sum()
        )

        correct = (
            (x["status"] == "C")
            & (x["pred_taxid"] == true_taxid)
        )

        n_correct = int(correct.sum())
        n_incorrect = n_reads - n_correct

        pct_correct = (
            100.0 * n_correct / n_reads
            if n_reads > 0
            else float("nan")
        )

        summary_rows.append(
            {
                "accession": accession,
                "dataset": dataset,
                "db": db,
                "true_taxid": true_taxid,
                "n_reads": n_reads,
                "n_classified": n_classified,
                "n_correct": n_correct,
                "n_incorrect": n_incorrect,
                "pct_correct": pct_correct,
                "file_found": True,
            }
        )

        # ----------------------------------------------------
        # Keep read ID + numeric predicted taxid
        # ----------------------------------------------------

        column_name = f"{dataset}_db{db}"

        x = x[
            ["read_id", "pred_taxid"]
        ].rename(
            columns={
                "pred_taxid": column_name
            }
        )

        # ----------------------------------------------------
        # Merge read classifications
        # ----------------------------------------------------

        if result is None:
            result = x

        else:
            old_n = len(result)

            result = result.merge(
                x,
                on="read_id",
                how="outer",
                validate="one_to_one",
                sort=False,
            )

            if len(result) != old_n:
                print(
                    f"WARNING: read sets differ after db{db}: "
                    f"{old_n:,} -> {len(result):,}"
                )

        del x

    # --------------------------------------------------------
    # If no files existed, stop with an error
    # --------------------------------------------------------

    if result is None:
        sys.exit(
            f"ERROR: no Kraken output files found for "
            f"{accession} ({dataset})"
        )

    # --------------------------------------------------------
    # Add true taxid
    # --------------------------------------------------------

    result.insert(
        1,
        "true_taxid",
        true_taxid,
    )

    expected_columns = [
        f"{dataset}_db{i}"
        for i in range(1, args.n_db + 1)
    ]

    present_columns = [
        col
        for col in expected_columns
        if col in result.columns
    ]

    result = result[
        ["read_id", "true_taxid"]
        + present_columns
    ]

    # --------------------------------------------------------
    # Write accession classification table
    # --------------------------------------------------------

    outfile = (
        out_dir
        / f"{accession}_{dataset}_classification.tsv"
    )

    result.to_csv(
        outfile,
        sep="\t",
        index=False,
    )

    print(f"Wrote: {outfile}")
    print(f"Reads: {len(result):,}")

    # --------------------------------------------------------
    # Write small summary fragment for this array task
    # --------------------------------------------------------

    summary = pd.DataFrame(summary_rows)

    summary_file = (
        summary_dir
        / f"{accession}_{dataset}_summary.tsv"
    )

    summary.to_csv(
        summary_file,
        sep="\t",
        index=False,
    )

    print(f"Wrote: {summary_file}")
    print("Done.")


if __name__ == "__main__":
    main()
