#!/usr/bin/env python3
# ABOUTME: Classifies detected mutant peptides as TP/FP against a cohort's known mutations.
# ABOUTME: Applies canonical-peptide filtering, then reports counts, thresholds, and per-sample breakdown.

# Cohort-specific behaviour comes from the optional TSVs run_samples.py reads, not from code:
# which injections share a sample, which truth-table names are spelled differently in the run
# names, and which samples are pools of others. With none of them the script does the generic
# thing -- one sample per run, names as written, no pools -- which is what an external user with
# their own data gets.

import argparse
import os
import sys
import re
import pyarrow.parquet as pq
import pandas as pd

import fragment_geometry as fg
import run_samples


# --- Configuration ---

# Detection gates (name, {column: max_value}), applied per precursor x run.
# LENIENT is the study-wide identification FDR, already enforced by DIA-NN on output.
# STRICT gates on Q.Value, which is recomputed inside each run. Lib.Q.Value carries one
# value per peptide for the whole study, so it cannot say whether the peptide was present
# in THIS sample -- and a replicate count built on it can only come out all-or-nothing.
# Built from the configured gate rather than fixed, so this script and gate_variant_cells.py
# cannot be handed different numbers and disagree about what was detected.
LENIENT_LIB_QVALUE = 0.01
QVALUE_GATE = 0.001


def thresholds(qvalue):
    return [
        ("UNFILTERED", {}),
        ("LENIENT", {"Lib.Q.Value": LENIENT_LIB_QVALUE}),
        ("STRICT",  {"Q.Value": qvalue}),
    ]

# Injections of a sample a call must be seen in. 1 means the requirement is off, which is the
# only honest default when the cohort's replicate structure is not described: see run_samples.py.
# The published cohorts used 2.
MIN_REPLICATES = 1

# Fraction of a variant peptide's library fragment intensity that must come from ions containing
# the mutated residue. A fragment that misses the site has the same mass in the mutant as in the
# wild type, so it is signal the normal protein produces identically and cannot be evidence for
# the mutation. Imported rather than declared so gate_variant_cells.py gates on the same value.
FRAGMENT_MIN = fg.FRAGMENT_SPECIFICITY_MIN

MUT_PATTERN = re.compile(r'^(.+?)_([A-Z*]\d+)_([A-Za-z*]+):(\d+)$')

# Indels carry the affected region and an HGVS change instead of a single ref/alt pair:
# P35222_23-71_D32_S33del:1, P04637_173-187_S183_G187delinsC:1, P51532_546_K546del:5.
# Without this the annotation does not parse, the detection gets no Gene, and groupby drops the
# row -- so an indel could never appear in the results table however good its evidence was.
# Tried after MUT_PATTERN, which is the stricter of the two.
INDEL_PATTERN = re.compile(r'^(.+?)_(\d+(?:-\d+)?)_(.+):(\d+)$')


def load_noncanonical_sequences(path):
    """Parse non_canonical_peptide_headers.txt → set of stripped sequences.

    Header format: PROTEINID_MUTATION:COUNT_SEQUENCE_CHARGE
    We extract the sequence (second-to-last _-delimited field).
    """
    seqs = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.rsplit('_', 2)
            if len(parts) >= 2:
                seqs.add(parts[-2])
    return seqs


def load_accession_gene_map(fasta_path):
    """Map UniProt accession -> gene symbol from the search FASTA's GN= field.

    Protein.Ids carries accessions; the ground-truth table carries gene symbols.
    """
    acc2gene = {}
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            fields = line[1:].split('|')
            if len(fields) < 2:
                continue
            m = re.search(r'GN=([^\s]+)', line)
            if m:
                # variant entries carry GN=KRAS_G12_V:657; keep the bare symbol
                acc2gene[fields[1]] = m.group(1).split('_')[0]
    return acc2gene


def parse_all_mutations(protein_ids, acc2gene):
    """Parse every mutation annotation in a semicolon-separated Protein.Ids field.

    Protein.Ids lists ALL proteins a peptide is consistent with, most recurrent first.
    The Genes field names only ONE arbitrary paralog, so a peptide shared between
    KRAS/NRAS/HRAS is labelled with whichever DIA-NN happened to pick -- which made
    genuine KRAS G12V and G12S calls score as false positives. Hence Protein.Ids.

    Substitutions and stop-gains carry a ref/alt pair; indels and delins carry an HGVS change
    that is used as written. No detected variant peptide mixes the two shapes in one Protein.Ids,
    so adding the second pattern cannot change a label the first one already produced.

    Returns list of (gene, protein_change, recurrence_count) tuples.
    """
    results = []
    for part in protein_ids.split(';'):
        part = part.strip()
        m = MUT_PATTERN.match(part)
        if m:
            change = f"p.{m.group(2)}{m.group(3).replace('*', 'Ter')}"
        else:
            m = INDEL_PATTERN.match(part)
            if not m:
                continue
            change = f"p.{m.group(3)}"
        acc = m.group(1)
        gene = acc2gene.get(acc) or acc2gene.get(acc.split('-')[0]) or acc
        results.append((gene, change, int(m.group(4))))
    return results


def parse_primary_mutation(muts, sample, ground_truth):
    """Reporting label: prefer a paralog the genome confirms, else the most recurrent.

    For a shared peptide every paralog is equally consistent with the spectrum, so the
    label is a naming choice, not a detection claim. Preferring the confirmed paralog
    keeps the output readable; the ambiguity is recorded in Paralog_Candidates.
    """
    if not muts:
        return None, None
    for gene, pc, _ in muts:
        if (sample, gene, pc) in ground_truth:
            return gene, pc
    best = max(muts, key=lambda x: x[2])
    return best[0], best[1]


def load_ground_truth(path, aliases):
    """Truth table → set of (sample, gene, protein_change) + the DataFrame.

    Sample names are put through the alias map on the way in, so every later lookup compares
    run-side spellings only. The DataFrame keeps a Sample_Matched column for the same reason,
    while Sample stays as written for reporting.
    """
    df = pd.read_csv(path, sep='\t')
    df['Sample_Matched'] = df['Sample'].map(lambda x: aliases.get(x, x))
    truth = set()
    for _, row in df.iterrows():
        truth.add((row['Sample_Matched'], row['Gene'], row['Protein.Change']))
    return truth, df


def classify(sample, all_mutations, ground_truth, all_truth_mutations, pools):
    """Classify a detection as TP, FP, POOL_TP or POOL_FP.

    A pool is a mixture, so no single sample's genome can confirm or refute a call in it: its
    hits are reported separately rather than counted as either. A pool with recorded members is
    checked against those members; one with '*' accepts any mutation the truth table knows,
    which is the conservative reading when the composition is not recorded.
    """
    if sample in pools:
        members = pools[sample]
        if members is None:
            for gene, pc in all_mutations:
                if (gene, pc) in all_truth_mutations:
                    return "POOL_TP"
            return "POOL_FP"
        for member in members:
            for gene, pc in all_mutations:
                if (member, gene, pc) in ground_truth:
                    return "POOL_TP"
        return "POOL_FP"

    for gene, pc in all_mutations:
        if (sample, gene, pc) in ground_truth:
            return "TP"
    return "FP"


def passes_threshold(row, thresh_dict):
    """Check if a row passes all conditions in a threshold dict."""
    return all(row[col] <= val for col, val in thresh_dict.items())


def main():
    parser = argparse.ArgumentParser(
        description="Classify detected mutant peptides as TP/FP against a cohort's known mutations"
    )
    script_dir = os.path.dirname(os.path.abspath(__file__))

    parser.add_argument(
        "parquet", nargs="?",
        default=os.path.join(script_dir, "Reports", "report_peptidoforms.parquet"),
        help="Path to report_peptidoforms.parquet",
    )
    parser.add_argument(
        "-t", "--truth", required=True,
        help="Ground-truth table: TSV with Sample, Gene, Protein.Change and "
             "Detected.By.DIANN columns. Required -- there is nothing to classify against "
             "without it. See data/truth/ for the two published cohorts",
    )
    parser.add_argument(
        "-s", "--sample-map", default="",
        help="TSV (run, sample, optional run_label) grouping injections into samples. "
             "Absent: every run is its own sample (default: none)",
    )
    parser.add_argument(
        "--aliases", default="",
        help="TSV (truth_name, run_name) for samples the truth table spells differently from "
             "the run names (default: none)",
    )
    parser.add_argument(
        "--pools", default="",
        help="TSV (pool, members) naming samples that are mixtures: '*' accepts any mutation "
             "in the truth table, empty accepts none, else 'A;B'. Hits in a pool are reported "
             "as POOL_TP/POOL_FP, never as TP/FP (default: none)",
    )
    parser.add_argument(
        "-q", "--qvalue", type=float, default=QVALUE_GATE,
        help="Per-run Q.Value a variant precursor must reach to count as present in that "
             "sample; the STRICT tier (default: %(default)s)",
    )
    parser.add_argument(
        "--fragment-min", type=float, default=FRAGMENT_MIN,
        help="Minimum share of library fragment intensity from ions containing the mutated "
             "residue (default: %(default)s)",
    )
    parser.add_argument(
        "--min-replicates", type=int, default=MIN_REPLICATES,
        help="Injections of a sample a call must pass the strict Q.Value gate in "
             "(default: %(default)s)",
    )
    parser.add_argument(
        "-n", "--non-canonical",
        default=os.path.join(script_dir, "non_canonical_peptide_headers.txt"),
        help="Path to non_canonical_peptide_headers.txt (canonical filter). "
             "Use --no-canonical-filter to skip.",
    )
    parser.add_argument(
        "--no-canonical-filter", action="store_true",
        help="Skip canonical peptide filtering (report all mutant precursors)",
    )
    parser.add_argument(
        "-p", "--proteome",
        default=os.path.join(script_dir, os.pardir, "data", "fasta", "proteome.fasta"),
        help="Search FASTA: maps UniProt accessions in Protein.Ids to gene symbols, and "
             "supplies the variant sequences the fragment filter locates the mutated residue in",
    )
    parser.add_argument(
        "--library",
        default=os.path.join(script_dir, "Library", "library_FROM_peptidoform.parquet"),
        help="DIA-NN library parquet, the source of the fragment intensities the "
             "fragment-geometry filter weighs (default: %(default)s)",
    )
    parser.add_argument(
        "-o", "--output",
        default=os.path.join(script_dir, "Peptidomics_Results",
                             "hotspot_detection_classification.tsv"),
        help="Output TSV path",
    )
    args = parser.parse_args()

    # ---- Load inputs ----
    aliases = run_samples.load_aliases(args.aliases)
    pools = run_samples.load_pools(args.pools)
    sample_map = run_samples.load_sample_map(args.sample_map)
    ground_truth, truth_df = load_ground_truth(args.truth, aliases)
    # A pool's own mutations are not knowable, so its hits are scored against every sample's.
    all_truth_mutations = {
        (row['Gene'], row['Protein.Change']) for _, row in truth_df.iterrows()
    }
    print(f"  truth table: {len(truth_df)} rows, {len(aliases)} alias(es), "
          f"{len(pools)} pool(s); replicate grouping "
          f"{'from ' + os.path.basename(args.sample_map) if sample_map else 'one sample per run'}",
          file=sys.stderr)

    read_cols = [
        "Run", "Precursor.Id", "Modified.Sequence", "Stripped.Sequence",
        "Precursor.Charge", "Protein.Ids", "Genes",
        "Precursor.Quantity",
        "Q.Value", "Lib.Q.Value", "Global.Q.Value", "PEP",
    ]
    print(f"Reading {args.parquet} ...", file=sys.stderr)
    df = pq.read_table(args.parquet, columns=read_cols).to_pandas()
    print(f"  Total rows: {len(df):,}", file=sys.stderr)

    # Filter for rows containing mutation annotations. Must test Protein.Ids, not Genes:
    # the Genes field omits the annotation for many rows that Protein.Ids carries.
    mut_mask = df["Protein.Ids"].str.contains(
        fg.VARIANT_ANNOTATION, na=False, regex=True
    )
    mut_df = df[mut_mask].copy()
    print(f"  Mutant precursor rows: {len(mut_df):,}", file=sys.stderr)
    n_genes_only = df["Genes"].str.contains(
        fg.VARIANT_ANNOTATION, na=False, regex=True
    ).sum()
    print(f"    (the Genes field would have found {n_genes_only:,})", file=sys.stderr)

    # Canonical peptide filtering
    if not args.no_canonical_filter:
        nc_seqs = load_noncanonical_sequences(args.non_canonical)
        n_before = len(mut_df)
        mut_df = mut_df[mut_df["Stripped.Sequence"].isin(nc_seqs)].copy()
        print(f"  Non-canonical sequences: {len(nc_seqs)}", file=sys.stderr)
        print(f"  After canonical filter: {len(mut_df):,} rows "
              f"(removed {n_before - len(mut_df):,} canonical matches)",
              file=sys.stderr)
    else:
        print("  Canonical filtering: SKIPPED", file=sys.stderr)

    # ---- Fragment geometry: can this peptide's fragments see its own mutation? ----
    # Scored per (sequence, charge) and never joined on a gene label: parse_primary_mutation
    # relabels a shared peptide to whichever paralog the genome confirms while the peptide list
    # keeps the accession DIA-NN assigned, so a (Gene, Protein.Change) join would silently void
    # every such call.
    frag_fractions, frag_skipped = fg.peptide_fractions(
        args.non_canonical, args.proteome, args.library
    )
    print(f"  Fragment specificity: scored {len(frag_fractions)} variant precursors "
          f"(gate >= {args.fragment_min:.0%})", file=sys.stderr)
    for header, reason in frag_skipped:
        print(f"    unscored: {header} -- {reason}", file=sys.stderr)
    mut_df["Frag_Specificity"] = [
        frag_fractions.get((seq, int(z)), float("nan"))
        for seq, z in zip(mut_df["Stripped.Sequence"], mut_df["Precursor.Charge"])
    ]

    # ---- Parse mutations, map to samples, classify ----
    acc2gene = load_accession_gene_map(args.proteome)
    print(f"  Accession->gene map: {len(acc2gene):,} entries", file=sys.stderr)

    mut_df["Sample"] = mut_df["Run"].apply(sample_map.sample)
    muts_per_row = [parse_all_mutations(p, acc2gene) for p in mut_df["Protein.Ids"]]

    # A row no pattern can read gets no Gene, and the groupby below then drops it from the
    # aggregation without trace. Name it here: an unreadable annotation must cost a visible
    # warning rather than silently shrink the results table.
    unlabelled = mut_df.loc[[not muts for muts in muts_per_row], "Protein.Ids"]
    if len(unlabelled):
        print(f"  WARNING: {len(unlabelled)} detection(s) carry an unreadable mutation "
              f"annotation and are excluded from the results table:", file=sys.stderr)
        for pid in sorted(set(unlabelled)):
            print(f"    {pid}", file=sys.stderr)

    labels = [
        parse_primary_mutation(muts, sample, ground_truth)
        for muts, sample in zip(muts_per_row, mut_df["Sample"])
    ]
    mut_df["Gene"] = [g for g, _ in labels]
    mut_df["Protein.Change"] = [pc for _, pc in labels]

    # Record the full ambiguity so a shared-peptide call is never read as gene-specific
    mut_df["Paralog_Candidates"] = [
        ";".join(sorted({f"{g}:{pc}" for g, pc, _ in muts})) for muts in muts_per_row
    ]
    mut_df["N_Paralog_Candidates"] = [
        len({g for g, _, _ in muts}) for muts in muts_per_row
    ]

    mut_df["Classification"] = [
        classify(sample, [(g, pc) for g, pc, _ in muts], ground_truth,
                 all_truth_mutations, pools)
        for sample, muts in zip(mut_df["Sample"], muts_per_row)
    ]

    # ---- Apply thresholds ----
    tiers = thresholds(args.qvalue)
    for name, thresh in tiers:
        mut_df[f"Pass_{name}"] = mut_df.apply(
            lambda row, t=thresh: passes_threshold(row, t), axis=1
        )

    # ---- Aggregate to mutation × sample level ----
    agg_records = []
    for (gene, pc, sample, cls), grp in mut_df.groupby(
        ["Gene", "Protein.Change", "Sample", "Classification"]
    ):
        # Reuse the per-row STRICT column so the replicate count and the STRICT tier
        # can never drift apart.
        strict_mask = grp["Pass_STRICT"]
        rec = {
            "Gene": gene,
            "Protein.Change": pc,
            "Sample": sample,
            "Classification": cls,
            "N_Paralog_Candidates": grp["N_Paralog_Candidates"].max(),
            "Paralog_Candidates": ";".join(sorted(set(
                c for s in grp["Paralog_Candidates"] for c in s.split(";") if c
            ))),
            "N_Precursors": grp["Precursor.Id"].nunique(),
            "N_Runs_Unfiltered": grp["Run"].nunique(),
            "N_Runs_STRICT": grp.loc[strict_mask, "Run"].nunique(),
            "N_Detections": len(grp),
            "Best_Q.Value": grp["Q.Value"].min(),
            "Best_Lib.Q.Value": grp["Lib.Q.Value"].min(),
            "Best_Global.Q.Value": grp["Global.Q.Value"].min(),
            "Best_PEP": grp.loc[strict_mask, "PEP"].min() if strict_mask.any() else float("nan"),
            # Max over the call's precursors: a call is as good as its best-placed peptide, and
            # charge states of the same peptide score differently.
            "Frag_Specificity": grp["Frag_Specificity"].max(),
        }
        for name, _ in tiers:
            col = f"Pass_{name}"
            rec[f"N_{col}"] = int(grp[col].sum())
            rec[f"Any_{col}"] = bool(grp[col].any())
        runs = sorted(grp["Run"].unique())
        rec["Runs"] = ";".join(sample_map.label(r) for r in runs)
        agg_records.append(rec)

    agg_df = pd.DataFrame(agg_records).sort_values(
        ["Gene", "Protein.Change", "Sample"]
    )

    # ---- Variant-specific filters (post-aggregation) ----
    # PEP is reported but not gated: it correlates 0.99 with Q.Value, so gating on both
    # filters the same evidence twice and costs true calls without removing false ones.
    agg_df["Pass_PEP_001"] = agg_df["Best_PEP"] <= 0.01
    agg_df["Pass_Replicate"] = agg_df["N_Runs_STRICT"] >= args.min_replicates
    # Frag_Specificity is reported whether it passes or not, so a rejection stays auditable.
    # A call with no score fails: an unparsable annotation is not evidence that the peptide can
    # carry its mutation. The count is printed below so this can never happen silently.
    agg_df["Pass_Fragment"] = agg_df["Frag_Specificity"] >= args.fragment_min
    agg_df["Pass_Variant_Filter"] = (
        agg_df["Any_Pass_STRICT"] & agg_df["Pass_Replicate"] & agg_df["Pass_Fragment"]
    )

    # ---- Write output TSV ----
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    agg_df.to_csv(args.output, sep="\t", index=False)
    print(f"\nDetailed table written to {args.output} ({len(agg_df)} rows)\n",
          file=sys.stderr)

    # ---- Print summary to stdout ----
    print_summary(mut_df, agg_df, truth_df, args.min_replicates, pools,
                  tiers, args.qvalue, args.fragment_min)


def print_summary(mut_df, agg_df, truth_df, min_replicates, pools,
                  tiers, qvalue, fragment_min):
    """Print human-readable classification report to stdout."""
    sep = "=" * 90
    # POOL_TP/POOL_FP are 4 characters longer than TP/FP, so the width is taken from the labels
    # actually present. With no pools this is 6 and the layout is the one without them.
    cls_w = max(6, *(len(c) for c in agg_df["Classification"].unique())) if len(agg_df) else 6

    # --- Section 1: Ground truth recall ---
    print(sep)
    print("SECTION 1: GROUND TRUTH RECALL")
    print(sep)
    print()
    print(f"{'Sample':<14s} {'Gene':<10s} {'Protein.Change':<14s} "
          f"{'Expected':<14s} {'Detected':<10s} {'Samples Detected In'}")
    print("-" * 90)

    for _, trow in truth_df.iterrows():
        sample = trow['Sample_Matched']
        gene = trow['Gene']
        pc = trow['Protein.Change']
        expected = "detectable" if trow['Detected.By.DIANN'] else "undetectable"

        match = agg_df[
            (agg_df['Gene'] == gene)
            & (agg_df['Protein.Change'] == pc)
            & (agg_df['Sample'] == sample)
        ]
        all_match = agg_df[
            (agg_df['Gene'] == gene) & (agg_df['Protein.Change'] == pc)
        ]
        if len(match) > 0:
            detected = "YES"
            samples_detected = ", ".join(sorted(all_match['Sample'].unique()))
        elif len(all_match) > 0:
            detected = "NO*"
            samples_detected = (
                f"(in other: {', '.join(sorted(all_match['Sample'].unique()))})"
            )
        else:
            detected = "NO"
            samples_detected = ""

        print(f"{trow['Sample']:<14s} {gene:<10s} {pc:<14s} "
              f"{expected:<14s} {detected:<10s} {samples_detected}")

    print()
    print("* NO* = not detected in the expected sample but detected in others")
    print()

    # --- Section 2: All detected mutations ---
    print(sep)
    print("SECTION 2: ALL POSITIVES DETECTED")
    print(sep)
    print()

    n_mutations = agg_df.groupby(["Gene", "Protein.Change"]).ngroups
    n_combos = len(agg_df)
    n_samples = agg_df["Sample"].nunique()

    print(f"Unique mutations detected:                  {n_mutations}")
    print(f"Unique mutation x sample combinations:       {n_combos}")
    print(f"Samples with detections:                     {n_samples}")
    print(f"Total run x precursor detections:            {len(mut_df)}")
    n_unlabelled = int(mut_df["Gene"].isna().sum())
    if n_unlabelled:
        print(f"Excluded, unreadable annotation:             {n_unlabelled}")
    print()

    print("-" * 90)
    print("Classification breakdown (mutation x sample level):")
    print("-" * 90)
    for cls in ["TP", "FP", "POOL_TP", "POOL_FP"]:
        subset = agg_df[agg_df["Classification"] == cls]
        if len(subset) == 0:
            continue
        n_m = subset.groupby(["Gene", "Protein.Change"]).ngroups
        n_s = subset["Sample"].nunique()
        print(f"  {cls:<10s}  {len(subset):>3d} combinations  "
              f"({n_m} unique mutations across {n_s} samples)")
    print()

    # --- Section 3: Threshold analysis ---
    print(sep)
    print("SECTION 3: THRESHOLD ANALYSIS")
    print(sep)
    print()
    print("Number of mutation x sample combinations with at least one "
          "precursor passing:")
    print()

    hdr = f"{'':14s}"
    for name, _ in tiers:
        hdr += f"  {name:>12s}"
    print(hdr)
    print("-" * (14 + 14 * len(tiers)))

    for cls in ["TP", "FP", "POOL_TP", "POOL_FP"]:
        subset = agg_df[agg_df["Classification"] == cls]
        if len(subset) == 0:
            continue
        line = f"  {cls:<12s}"
        for name, _ in tiers:
            n_pass = int(subset[f"Any_Pass_{name}"].sum())
            line += f"  {n_pass:>12d}"
        print(line)

    line = f"  {'TOTAL':<12s}"
    for name, _ in tiers:
        n_pass = int(agg_df[f"Any_Pass_{name}"].sum())
        line += f"  {n_pass:>12d}"
    print(line)
    print()

    # --- Section 3b: Variant-specific filter analysis ---
    print(sep)
    print("SECTION 3b: VARIANT-SPECIFIC FILTER ANALYSIS")
    print(sep)
    print()
    print("Filters applied at mutation x sample level (post-aggregation):")
    print("  Pass_PEP_001:        Best_PEP (from STRICT-passing precursors) <= 0.01 "
          "[reported, not gated]")
    print(f"  Pass_Replicate:      N_Runs_STRICT >= {min_replicates} "
          f"(Q.Value <= {qvalue} in that many injections of the sample)")
    print(f"  Pass_Fragment:       Frag_Specificity >= {fragment_min:.0%} (share of library "
          "fragment intensity from ions containing the mutated residue)")
    print("  Pass_Variant_Filter: Any_Pass_STRICT AND Pass_Replicate AND Pass_Fragment")
    print()
    filt_hdr = (f"{'':14s}  {'Pass_PEP_001':>12s}  {'Pass_Replicate':>14s}"
                f"  {'Pass_Fragment':>13s}  {'Pass_Variant_Filter':>19s}")
    print(filt_hdr)
    print("-" * 83)
    for cls in ["TP", "FP", "POOL_TP", "POOL_FP"]:
        subset = agg_df[agg_df["Classification"] == cls]
        if len(subset) == 0:
            continue
        n_pep = int(subset["Pass_PEP_001"].sum())
        n_rep = int(subset["Pass_Replicate"].sum())
        n_frag = int(subset["Pass_Fragment"].sum())
        n_both = int(subset["Pass_Variant_Filter"].sum())
        print(f"  {cls:<12s}  {n_pep:>12d}  {n_rep:>14d}  {n_frag:>13d}  {n_both:>19d}")
    n_pep = int(agg_df["Pass_PEP_001"].sum())
    n_rep = int(agg_df["Pass_Replicate"].sum())
    n_frag = int(agg_df["Pass_Fragment"].sum())
    n_both = int(agg_df["Pass_Variant_Filter"].sum())
    print(f"  {'TOTAL':<12s}  {n_pep:>12d}  {n_rep:>14d}  {n_frag:>13d}  {n_both:>19d}")
    n_unscored = int(agg_df["Frag_Specificity"].isna().sum())
    if n_unscored:
        print(f"\n  {n_unscored} call(s) carry no fragment score and so fail Pass_Fragment "
              "(see the 'unscored' lines on stderr)")
    print()

    # --- Section 4: Per-mutation detail ---
    print(sep)
    print("SECTION 4: PER-MUTATION DETAIL")
    print(sep)
    print()

    for (gene, pc), grp in agg_df.groupby(["Gene", "Protein.Change"]):
        truth_rows = truth_df[
            (truth_df['Gene'] == gene) & (truth_df['Protein.Change'] == pc)
        ]
        if len(truth_rows) > 0:
            gt_parts = []
            for _, tr in truth_rows.iterrows():
                tag = "detectable" if tr['Detected.By.DIANN'] else "undetectable"
                # The run-side spelling, so a reader can find the sample in the rows below.
                gt_parts.append(f"{tr['Sample_Matched']}({tag})")
            gt_str = f"in the truth table: {', '.join(gt_parts)}"
        else:
            gt_str = "NOT in the truth table"

        print(f"  {gene} {pc}  [{gt_str}]")
        for _, row in grp.iterrows():
            best_thresh = "NONE"
            for name, _ in reversed(tiers):
                if name == "UNFILTERED":
                    continue
                if row[f"Any_Pass_{name}"]:
                    best_thresh = name
                    break

            pep_str = f"{row['Best_PEP']:.4f}" if pd.notna(row['Best_PEP']) else "N/A"
            frag_str = (f"{row['Frag_Specificity']:.1%}"
                        if pd.notna(row['Frag_Specificity']) else "N/A")
            print(f"    {row['Sample']:<14s}  {row['Classification']:<{cls_w}s}  "
                  f"{row['N_Precursors']:>2d} precursors  "
                  f"{row['N_Runs_STRICT']:>2d}/{row['N_Runs_Unfiltered']:>2d} runs(strict/total)  "
                  f"bestQ={row['Best_Q.Value']:.4f}  "
                  f"bestLibQ={row['Best_Lib.Q.Value']:.4f}  "
                  f"PEP={pep_str}  "
                  f"frag={frag_str:>6s}  "
                  f"[{best_thresh}]")
        print()

    # --- Section 5: Per-sample summary ---
    print(sep)
    print("SECTION 5: PER-SAMPLE SUMMARY")
    print(sep)
    print()

    for sample in sorted(agg_df["Sample"].unique()):
        s_data = agg_df[agg_df["Sample"] == sample]
        n_strict = int(s_data["Any_Pass_STRICT"].sum())
        # A pool's hits are POOL_TP/POOL_FP and are counted as such. Folding them into TP/FP
        # here would contradict the rest of the report, which never calls a pool hit either:
        # no single sample's genome can confirm or refute a call in a mixture.
        if sample in pools:
            counted = (f"{int((s_data['Classification'] == 'POOL_TP').sum())} POOL_TP, "
                       f"{int((s_data['Classification'] == 'POOL_FP').sum())} POOL_FP")
        else:
            counted = (f"{int((s_data['Classification'] == 'TP').sum())} TP, "
                       f"{int((s_data['Classification'] == 'FP').sum())} FP")

        print(f"  {sample}: {len(s_data)} mutations detected  "
              f"({counted})  "
              f"[{n_strict} pass STRICT]")
        for _, row in s_data.sort_values("Classification").iterrows():
            print(f"    {row['Classification']:<{cls_w}s}  "
                  f"{row['Gene']:<10s} {row['Protein.Change']:<14s}  "
                  f"{row['N_Precursors']} precursors  "
                  f"bestQ={row['Best_Q.Value']:.4f}")
    print()


if __name__ == "__main__":
    main()
