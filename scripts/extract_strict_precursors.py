#!/usr/bin/env python3
# ABOUTME: Extracts precursors from DIA-NN peptidoform report at strict q-value thresholds.
# ABOUTME: Designed for variant peptide detection pipeline — outputs filtered precursors for downstream canonical filtering.

import argparse
import os
import sys
import pyarrow.parquet as pq
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Extract precursors passing strict q-value thresholds from DIA-NN peptidoform report"
    )
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_input = os.path.join(script_dir, "Reports", "report_peptidoforms.parquet")
    default_output = os.path.join(script_dir, "Reports", "strict_precursors_peptide_list.tsv")

    parser.add_argument(
        "input", nargs="?", default=default_input,
        help="Path to report_peptidoforms.parquet (default: %(default)s)"
    )
    parser.add_argument(
        "-o", "--output", default=default_output,
        help="Output TSV path (default: %(default)s)"
    )
    parser.add_argument(
        "--lib-qvalue", type=float, default=0.001,
        help="Lib.Q.Value threshold (default: 0.001 = 0.1%% FDR). "
             "With MBR enabled, Lib.* columns are statistically correct per DIA-NN guidance."
    )
    parser.add_argument(
        "--full-table", action="store_true", default=False,
        help="Also write the full per-run table (large; off by default)"
    )
    args = parser.parse_args()

    cols_output = [
        "Run",
        "Precursor.Id",
        "Modified.Sequence",
        "Stripped.Sequence",
        "Precursor.Charge",
        "Precursor.Mz",
        "Protein.Ids",
        "Protein.Group",
        "Protein.Names",
        "Genes",
        "Precursor.Quantity",
        "Precursor.Normalised",
        "RT",
        "Q.Value",
        "Lib.Q.Value",
        "Global.Q.Value",
        "PG.Q.Value",
        "GG.Q.Value",
    ]

    print(f"Reading {args.input} ...", file=sys.stderr)
    df = pq.read_table(args.input, columns=cols_output).to_pandas()
    n_total = len(df)
    print(f"  Total rows: {n_total:,}", file=sys.stderr)

    if n_total == 0:
        print("WARNING: No rows in input. DIA-NN may have failed to process files.", file=sys.stderr)
        pd.DataFrame(columns=["Stripped.Sequence"]).to_csv(args.output, sep="\t", index=False)
        sys.exit(0)

    if (df["Lib.Q.Value"] == 0).all():
        print("WARNING: All Lib.Q.Value entries are 0. This typically occurs in single-file "
              "runs where DIA-NN disables MBR. The strict filter will pass all rows. "
              "For meaningful FDR filtering, run with multiple input files.", file=sys.stderr)

    # --- Filtering ---
    # Lib.Q.Value is the statistically correct filter with MBR (--reanalyse).
    # Per DIA-NN guidance, Lib.* columns account for match-between-runs correctly.
    lq = args.lib_qvalue
    mask = df["Lib.Q.Value"] <= lq

    df_filt = df.loc[mask].copy()
    n_pass = len(df_filt)
    print(f"  Threshold: Lib.Q.Value <= {lq}",
          file=sys.stderr)
    print(f"  Rows passing: {n_pass:,} ({n_pass / n_total * 100:.1f}%)", file=sys.stderr)

    # Summary stats
    unique_prec = df_filt["Precursor.Id"].nunique()
    unique_pep = df_filt["Stripped.Sequence"].nunique()
    unique_proteins = df_filt["Protein.Ids"].nunique()
    n_runs = df_filt["Run"].nunique()
    print(f"  Unique precursors: {unique_prec:,}", file=sys.stderr)
    print(f"  Unique stripped peptides: {unique_pep:,}", file=sys.stderr)
    print(f"  Unique protein entries: {unique_proteins:,}", file=sys.stderr)
    print(f"  Runs represented: {n_runs:,}", file=sys.stderr)

    # Deduplicated peptide list (one row per unique stripped sequence)
    df_pep = (
        df_filt
        .groupby("Stripped.Sequence", sort=False)
        .agg(
            Modified_Sequences=("Modified.Sequence", lambda x: ";".join(sorted(x.unique()))),
            Protein_Ids=("Protein.Ids", lambda x: ";".join(sorted(x.unique()))),
            Genes=("Genes", lambda x: ";".join(sorted(x.unique()))),
            Runs_Detected=("Run", "nunique"),
            Best_Q_Value=("Q.Value", "min"),
            Best_Lib_Q_Value=("Lib.Q.Value", "min"),
        )
        .reset_index()
        .sort_values("Stripped.Sequence")
    )
    df_pep.to_csv(args.output, sep="\t", index=False)
    print(f"  Peptide list written to {args.output} ({len(df_pep):,} unique peptides)", file=sys.stderr)

    # Optional: full per-run table
    if args.full_table:
        full_path = args.output.replace("_peptide_list.tsv", ".tsv")
        if full_path == args.output:
            full_path = args.output.replace(".tsv", "_full.tsv")
        df_filt.sort_values(
            ["Run", "Protein.Ids", "Stripped.Sequence", "Precursor.Charge"],
            inplace=True,
        )
        df_filt.to_csv(full_path, sep="\t", index=False)
        print(f"  Full table written to {full_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
