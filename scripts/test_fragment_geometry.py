#!/usr/bin/env python3
# ABOUTME: Self-contained checks for the fragment-geometry rule -- no search output needed.
# ABOUTME: Run with: python3 scripts/test_fragment_geometry.py

# The rule derives the diagnostic fragment set from where a variant sequence stops matching its
# wild-type parent, which is what lets substitutions, insertions, deletions, delins and stop-gains
# score by the same arithmetic. These checks pin that derivation for one example of each class,
# using sequences written out here rather than a cohort's FASTA and library, so they run on a
# fresh clone before any data has been downloaded.
#
# The full regression -- every published substitution score held to nine decimal places, and the
# two indel scores -- needs a completed run's library_FROM_peptidoform.parquet and lives with the
# analysis, not in this repo. What is here is the part that can be checked from the rule alone.

import sys

import fragment_geometry as fg


def test_substitution_reduces_to_the_older_form():
    """For a substitution at peptide position p, (L, S) must equal (p-1, n-p).

    That is the b_k >= p / y_k >= n-p+1 rule the sequence-divergence form replaced, so any
    substitution scores exactly what it scored before.
    """
    wt = "MAAAKGGGGRSSSSK"
    variant = wt[:9] + "Q" + wt[10:]      # R10Q, so the peptide loses a cleavage site
    peptide = "GGGGQSSSSK"                # p = 5 within the peptide, n = 10
    assert variant.find(peptide) == 5, variant
    L, S = fg.wt_identical_flanks(variant, wt, peptide)
    p, n = 5, len(peptide)
    assert (L, S) == (p - 1, n - p), (L, S)
    return 1


def test_insertion_deletion_and_delins():
    """The junction has to fall strictly inside the peptide, whatever its length change."""
    # Deletion of PP: the divergent region is empty, so L + S == n and only the ions crossing
    # the junction differ.
    L, S = fg.wt_identical_flanks("MKLQWERTY", "MKLPPQWERTY", "LQWE")
    assert (L, S) == (1, 3), (L, S)
    # Insertion of PP: same junction test in the other direction.
    L, S = fg.wt_identical_flanks("MKLPPQWERTY", "MKLQWERTY", "LPPQWE")
    assert L == 1 and S == 3, (L, S)
    # Delins: two residues replaced by one. Here the divergent region is one residue wide, so
    # L + S == n - 1, unlike the pure deletion above where it is empty and L + S == n.
    L, S = fg.wt_identical_flanks("MKLCWERTY", "MKLQQWERTY", "LCWE")
    assert (L, S) == (1, 2), (L, S)
    return 3


def test_peptide_clear_of_the_variant_scores_nothing():
    """A peptide not reaching the divergence is wild type, on either side of it."""
    variant, wt = "AAAAACAAAAA", "AAAAAAAAAAA"
    assert fg.wt_identical_flanks(variant, wt, variant) == (5, 5)
    assert fg.wt_identical_flanks(variant, wt, "NOTPRESENT") is None
    # A peptide from before the substitution, in a sequence where it is locatable.
    assert fg.wt_identical_flanks("MKGGGGKCSSSS", "MKGGGGKASSSS", "GGGGK") is None
    return 3


def test_stop_gain_scores_on_its_new_c_terminus():
    """A truncated protein has no divergent residue: every y-ion is diagnostic, no b-ion is."""
    wt = "MKAAAAGGGGSSSSDDDD"
    variant = "MKAAAAGGGGS"        # stop after S, which is not K/R
    peptide = "AAAAGGGGS"          # ends at the new C-terminus
    assert fg.wt_identical_flanks(variant, wt, peptide) == (len(peptide), 0)
    # A peptide upstream of the truncation is pure wild-type sequence.
    assert fg.wt_identical_flanks(variant, wt, "MKAAAA") is None
    return 2


def test_stop_gain_at_a_cleavage_site_is_not_evidence():
    """A truncation landing straight after K/R leaves the wild type making the same peptide.

    Trypsin cuts there in both proteins, so not even the C-terminus is novel.
    """
    wt = "MKAAAAGGGGKSSSSDDDD"
    variant = "MKAAAAGGGGK"        # stop after K -- a site the wild type cleaves anyway
    assert fg.wt_identical_flanks(variant, wt, "AAAAGGGGK") is None
    return 1


def test_isoform_keeps_its_own_accession():
    """A variant built on an isoform must be compared against that isoform, not the canonical."""
    assert fg.base_accession("Q01196-8_R201_Q:2") == "Q01196-8"
    assert fg.base_accession("P35222_23-71_D32_S33del:1") == "P35222"
    assert fg.base_accession("P60484_R130_*:38") == "P60484"
    return 3


def test_spanning_fraction_weighs_by_intensity():
    """b_k counts when k > L, y_k when k > S; the denominator is the precursor's total."""
    import pandas as pd
    frags = pd.DataFrame({
        "Fragment.Type": ["b", "b", "y", "y"],
        "Fragment.Series.Number": [2, 5, 2, 5],
        "Relative.Intensity": [0.1, 0.2, 0.3, 0.4],
    })
    # L = 3, S = 3: b5 and y5 are diagnostic, b2 and y2 are not.
    got = fg.spanning_fraction(frags, 3, 3)
    assert abs(got - 0.6) < 1e-12, got
    # Nothing diagnostic at all.
    assert fg.spanning_fraction(frags, 9, 9) == 0.0
    # A precursor with no intensity cannot be scored, and must not divide by zero.
    frags["Relative.Intensity"] = 0.0
    assert fg.spanning_fraction(frags, 3, 3) is None
    return 3


def test_threshold_selects_the_failures():
    """failing_sequences returns the peptides below the cut, keyed on sequence alone."""
    fractions = {("PEPTIDEA", 2): 0.99, ("PEPTIDEB", 2): 0.05,
                 ("PEPTIDEC", 3): 0.15, ("PEPTIDEC", 2): 0.14}
    below = fg.failing_sequences(fractions, 0.15)
    assert "PEPTIDEB" in below and "PEPTIDEA" not in below, below
    # PEPTIDEC passes at one charge and fails at another; a charge state is its own library entry,
    # so the passing one keeps the peptide.
    assert "PEPTIDEC" not in below, below
    return 3


def test_variant_annotation_matches_every_class():
    """A report row the row-selection pattern misses is dropped before every later filter."""
    import re
    pat = re.compile(fg.VARIANT_ANNOTATION)
    for acc in ["P01111_Q61_R:204", "P60484_R130_*:38",
                "P35222_23-71_D32_S33del:1", "P04637_173-187_S183_G187delinsC:1",
                "P51532_546_K546del:5"]:
        assert pat.search(acc), acc
    for acc in ["P01111", "sp|P04637|P53_HUMAN", "Q01196-8"]:
        assert not pat.search(acc), acc
    return 8


if __name__ == "__main__":
    failed = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith("test_"):
            continue
        try:
            n = fn()
            print(f"PASS  {name}  ({n} checks)")
        except (AssertionError, SystemExit) as exc:
            # SystemExit too: the loaders exit on a malformed file, and an unexpected one must be
            # reported as this test failing rather than silently ending the whole run.
            failed += 1
            print(f"FAIL  {name}\n      {exc}")
    print("\nall checks passed" if not failed else f"\n{failed} test(s) failed")
    sys.exit(1 if failed else 0)
