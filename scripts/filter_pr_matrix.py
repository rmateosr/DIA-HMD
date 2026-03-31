#!/usr/bin/env python3
# ABOUTME: Filters a DIA-NN pr_matrix to retain only precursors whose stripped sequence
# ABOUTME: appears in a strict-filtered peptide list.

import argparse
import os
import sys
import pandas as pd


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))

    default_matrix = os.path.join(script_dir, "Reports", "report_peptidoforms.pr_matrix.tsv")
    default_peptides = os.path.join(script_dir, "Reports", "strict_precursors_peptide_list.tsv")
    default_output = os.path.join(script_dir, "Reports", "report_peptidoforms.pr_matrix.strict.tsv")

    parser = argparse.ArgumentParser(
        description="Filter pr_matrix to rows matching a strict-filtered peptide list"
    )
    parser.add_argument(
        "-m", "--matrix", default=default_matrix,
        help="Path to report_peptidoforms.pr_matrix.tsv (default: %(default)s)"
    )
    parser.add_argument(
        "-p", "--peptides", default=default_peptides,
        help="Path to strict_precursors_peptide_list.tsv (default: %(default)s)"
    )
    parser.add_argument(
        "-o", "--output", default=default_output,
        help="Output TSV path (default: %(default)s)"
    )
    args = parser.parse_args()

    print(f"Reading peptide list {args.peptides} ...", file=sys.stderr)
    keep = set(pd.read_csv(args.peptides, sep="\t", usecols=["Stripped.Sequence"])["Stripped.Sequence"])
    print(f"  {len(keep):,} unique peptides to keep", file=sys.stderr)

    print(f"Reading matrix {args.matrix} ...", file=sys.stderr)
    df = pd.read_csv(args.matrix, sep="\t", low_memory=False)
    n_total = len(df)
    print(f"  Total rows: {n_total:,}", file=sys.stderr)

    mask = df["Stripped.Sequence"].isin(keep)
    df_filt = df.loc[mask]
    n_pass = len(df_filt)
    print(f"  Rows passing filter: {n_pass:,} ({n_pass / n_total * 100:.1f}%)", file=sys.stderr)

    df_filt.to_csv(args.output, sep="\t", index=False)
    print(f"  Filtered matrix written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
