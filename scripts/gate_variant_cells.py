#!/usr/bin/env python3
# ABOUTME: Blanks per-run cells of variant-peptide rows in a DIA-NN pr_matrix that the classifier
# ABOUTME: would not report -- failing the per-run Q.Value, the replicate count, or fragment geometry.

# The wide pr_matrix holds intensities but no q-values, so this filter cannot be applied
# anywhere downstream of it -- the q-values exist only in the parquet report.
# Only variant-peptide rows are filtered: the wild-type counterpart rows are reference
# context, not detection claims, and their weaker charge states legitimately fail this
# threshold. Cells are emptied rather than rows deleted, so every metadata column survives
# and peptide.fasta and the list of peptides considered are unaffected.
#
# All three gates mirror classify_hotspot_detections.py, so a figure cannot plot a point the
# results table rejected. They differ in grain:
#   Q.Value    per cell     -- a peptide can be present in one sample and absent in another
#   replicate  per sample   -- a detection in fewer injections than --min-replicates is not
#                             reported, so the surviving cells are emptied too; otherwise a point
#                             is drawn for a call the results table rejected. Off by default
#   fragment   per peptide  -- a peptide whose library fragments cannot tell the mutant from the
#                             wild type is not evidence in any run, so every cell goes
# A variant peptide with no fragment score at all (an indel annotation the rule cannot parse) is
# emptied for the same reason: it produces no call in the results table either.
#
# The replicate gate is applied per precursor, whereas the classifier counts runs per call
# aggregated over all of a call's precursors. The two agree whenever some single precursor already
# reaches the required number of runs, which holds for every reported call in both cohorts. Where
# they could differ the script says so rather than diverging silently.

import argparse
import os
import sys

import pandas as pd
import pyarrow.parquet as pq

import fragment_geometry as fg
import run_samples

RUN_COL_MARKER = "raw.dia"

# Injections of a sample a variant precursor must pass the Q.Value gate in. 1 means the gate is
# off: with no sample map the cohort's replicate structure is unknown, and rejecting a peptide for
# appearing once would then reject every peptide. The published cohorts used 2 -- see the README.
MIN_REPLICATES = 1


def load_variant_sequences(path):
    """non_canonical_peptide_headers.txt -> set of stripped sequences.

    Header format: PROTEINID_MUTATION:COUNT_SEQUENCE_CHARGE, so the sequence is the
    second-to-last _-delimited field (the protein ID may itself contain underscores).
    """
    seqs = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.rsplit("_", 2)
            if len(parts) >= 2:
                seqs.add(parts[-2])
    return seqs


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    reports = os.path.join(script_dir, "Reports")

    parser = argparse.ArgumentParser(
        description="Blank pr_matrix cells whose per-run Q.Value fails the detection threshold"
    )
    parser.add_argument(
        "-m", "--matrix",
        default=os.path.join(reports, "report_peptidoforms.pr_matrix.tsv"),
        help="DIA-NN's own precursor matrix, at its 1%% FDR output (default: %(default)s)",
    )
    parser.add_argument(
        "-r", "--report", default=os.path.join(reports, "report_peptidoforms.parquet"),
        help="DIA-NN parquet report, the only source of per-run Q.Value (default: %(default)s)",
    )
    parser.add_argument(
        "-n", "--variant-peptides",
        default=os.path.join(script_dir, "non_canonical_peptide_headers.txt"),
        help="Variant peptide list; only these rows are gated (default: %(default)s)",
    )
    parser.add_argument(
        "-o", "--output",
        default=os.path.join(reports, "report_peptidoforms.pr_matrix.strict.tsv"),
        help="Output matrix (default: %(default)s)",
    )
    parser.add_argument(
        "-q", "--qvalue", type=float, default=0.001,
        help="Per-run Q.Value threshold for calling a peptide present in a sample "
             "(default: %(default)s)",
    )
    parser.add_argument(
        "-p", "--proteome",
        default=os.path.join(script_dir, os.pardir, "data", "fasta", "proteome.fasta"),
        help="Search FASTA, supplying the variant sequences the fragment filter locates the "
             "mutated residue in (default: %(default)s)",
    )
    parser.add_argument(
        "--library",
        default=os.path.join(script_dir, "Library", "library_FROM_peptidoform.parquet"),
        help="DIA-NN library parquet, the source of the fragment intensities the "
             "fragment-geometry filter weighs (default: %(default)s)",
    )
    parser.add_argument(
        "--fragment-min", type=float, default=fg.FRAGMENT_SPECIFICITY_MIN,
        help="Minimum share of library fragment intensity from ions containing the mutated "
             "residue (default: %(default)s)",
    )
    parser.add_argument(
        "--min-replicates", type=int, default=MIN_REPLICATES,
        help="Injections of a sample a variant precursor must pass the Q.Value gate in "
             "(default: %(default)s)",
    )
    parser.add_argument(
        "-s", "--sample-map", default="",
        help="TSV (run, sample, optional run_label) grouping injections into samples so the "
             "replicate requirement has something to count. Absent: every run is its own "
             "sample (default: none)",
    )
    args = parser.parse_args()

    variant_seqs = load_variant_sequences(args.variant_peptides)
    print(f"Variant peptides to gate: {len(variant_seqs)}", file=sys.stderr)

    frag_fractions, frag_skipped = fg.peptide_fractions(
        args.variant_peptides, args.proteome, args.library
    )
    for header, reason in frag_skipped:
        print(f"  unscored: {header} -- {reason}", file=sys.stderr)
    below = fg.failing_sequences(frag_fractions, args.fragment_min)
    unscored = variant_seqs - {peptide for peptide, _charge in frag_fractions}
    frag_rejected = below | unscored
    print(f"  fragment gate >= {args.fragment_min:.0%}: {len(below)} peptides below it, "
          f"{len(unscored)} unscored, {len(frag_rejected)} to be emptied in every run",
          file=sys.stderr)

    print(f"Reading {args.matrix} ...", file=sys.stderr)
    mx = pd.read_csv(args.matrix, sep="\t", low_memory=False)
    run_cols = [c for c in mx.columns if RUN_COL_MARKER in c]
    var_rows = mx.index[mx["Stripped.Sequence"].isin(variant_seqs)]
    frag_fail_rows = set(var_rows[mx.loc[var_rows, "Stripped.Sequence"].isin(frag_rejected)])
    print(f"  {len(mx):,} rows, {len(run_cols)} runs, {len(var_rows)} variant rows "
          f"({len(frag_fail_rows)} rejected by fragment geometry)", file=sys.stderr)
    if not len(var_rows):
        print("  WARNING: no variant rows matched -- check the peptide list", file=sys.stderr)

    print(f"Reading {args.report} ...", file=sys.stderr)
    rep = pq.read_table(args.report,
                        columns=["Run", "Precursor.Id", "Q.Value"]).to_pandas()
    var_ids = set(mx.loc[var_rows, "Precursor.Id"])
    rep = rep[rep["Precursor.Id"].isin(var_ids)]
    # Keyed on the bare run stem, because the parquet Run and the matrix column spell the same
    # run differently ('X.raw' against '/path/X.raw.dia').
    passing = set(
        zip(rep.loc[rep["Q.Value"] <= args.qvalue, "Run"].map(run_samples.run_key),
            rep.loc[rep["Q.Value"] <= args.qvalue, "Precursor.Id"])
    )
    print(f"  variant precursor x run detections passing Q <= {args.qvalue}: {len(passing):,}",
          file=sys.stderr)

    sample_map = run_samples.load_sample_map(args.sample_map)
    print(f"  replicate grouping: "
          f"{'sample map ' + os.path.basename(args.sample_map) if sample_map else 'one sample per run'}"
          f", requiring {args.min_replicates} injection(s)", file=sys.stderr)
    by_sample = {}
    for col in run_cols:
        by_sample.setdefault(sample_map.sample(col), []).append(col)

    blanked_q = blanked_frag = blanked_rep = kept = 0
    for i in var_rows:
        pid = mx.at[i, "Precursor.Id"]
        for cols_of_sample in by_sample.values():
            present = [c for c in cols_of_sample if pd.notna(mx.at[i, c])]
            if not present:
                continue
            # Fragment geometry first, so a cell failing several tests is attributed to the
            # broadest reason rather than counted twice.
            if i in frag_fail_rows:
                for c in present:
                    mx.at[i, c] = pd.NA
                blanked_frag += len(present)
                continue
            survives_q = [c for c in present if (run_samples.run_key(c), pid) in passing]
            for c in present:
                if c not in survives_q:
                    mx.at[i, c] = pd.NA
            blanked_q += len(present) - len(survives_q)
            # A detection in fewer injections than required is not reported, so its surviving
            # cells are emptied too. At the default of 1 this can never fire.
            if len(survives_q) < args.min_replicates:
                for c in survives_q:
                    mx.at[i, c] = pd.NA
                blanked_rep += len(survives_q)
            else:
                kept += len(survives_q)

    blanked = blanked_q + blanked_rep + blanked_frag
    total = blanked + kept
    print(f"  variant cells: {kept:,} kept, {blanked:,} blanked "
          f"({100 * blanked / total if total else 0:.1f}% of {total:,}) -- "
          f"{blanked_q:,} on Q.Value, {blanked_rep:,} on the "
          f"{args.min_replicates}-injection requirement, "
          f"{blanked_frag:,} on fragment geometry", file=sys.stderr)

    tmp = args.output + ".tmp"
    mx.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, args.output)
    print(f"  gated matrix written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
