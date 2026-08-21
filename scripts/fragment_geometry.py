#!/usr/bin/env python3
# ABOUTME: Fragment-geometry rule -- which of a variant peptide's library fragments differ in mass
# ABOUTME: from the wild type, and what fraction of the library fragment intensity they carry.

# A fragment whose mass is the same in the mutant as in the wild type cannot be evidence that the
# variant is present: it is signal the normal protein produces identically. Which fragments those
# are follows from where the variant sequence stops matching its wild-type parent. With the two
# proteins sharing a prefix of E residues and a suffix of F residues, and the peptide sitting at
# offset o with length n,
#     L = clamp(E - o, 0, n)                          residues still wild type at the N-terminal
#     S = clamp((o + n) - (len(variant) - F), 0, n)   and C-terminal end of the peptide
# and b_k differs from wild type iff k > L, y_k iff k > S. Ion type and series number are all this
# needs -- no residue masses, modifications or isotopes.
#
# Taking L and S from the sequences rather than from the accession's mutation annotation is what
# lets substitutions, insertions, deletions, delins and stop-gains score by the same arithmetic.
# For a substitution at peptide position p it reduces to L = p-1, S = n-p, i.e. the older
# b_k >= p / y_k >= n-p+1 form, and reproduces every score that form produced.
#
# The wild-type sequence must come from the search FASTA rather than a canonical proteome. A
# variant built on an isoform (Q01196-8) compared against the canonical sequence reads every
# isoform difference as part of the mutation; that alone moved ITVDGPQEPR from 0.967 to 0.985.
#
# The metric is a property of the library entry, so it is one value per precursor for the whole
# study, the same structural limitation as Lib.Q.Value. It can remove peptides incapable of
# proving their variant; it can never say why a peptide is right in sample A and wrong in B.
#
# Charge states are separate library rows and are counted separately, so the denominator is a sum
# of ratios rather than a physical intensity (y7 at 1+ and 2+ can together carry ~46% of one
# peptide's total). Collapsing them to one entry per (type, series) would change every value and
# has not been tested.
#
# A mass-arithmetic form of the rule was evaluated and rejected: it needs a residue-mass table
# plus fixed-modification handling and changes no result, because its only extra catch is I<->L,
# which the canonical filter already removes.

import re

import pyarrow.parquet as pq

# Minimum fraction of library fragment intensity that must come from ions whose mass the variant
# changes. 15% is the largest cost-free cut available in both cohorts: it drops no genomically
# confirmed call while taking cell-line precision 48.0% -> 80.0% and PDX 66.7% -> 81.5%. What caps
# it is the cleavage-removing class (ref = K/R), whose mutant peptide fuses two wild-type peptides
# so only the weak long b-ions are diagnostic -- its lowest true call, PIK3R1 p.K567E, sits at
# 17.9%. The value was fitted on substitutions only; indels are scored against it unchanged.
FRAGMENT_SPECIFICITY_MIN = 0.15

# Variant accession, e.g. P04637_R181_H:9 -> base P04637, ref R, protein position 181, alt H.
# Substitutions and stop-gains only. Kept because pool_dilution.py matches against it directly;
# the scoring below no longer needs it.
VARIANT_ACCESSION = re.compile(r'^(.+?)_([A-Z*]\d+)_([A-Za-z*]+):(\d+)$')

# Row-selection form of the same idea, for scanning a semicolon-separated Protein.Ids field:
# a substitution or stop-gain (P01111_Q61_R:204) or an indel (P35222_23-71_D32_S33del:1).
# Shared so the classifier and the cascade report cannot disagree about which report rows carry a
# variant at all. A row this misses is dropped before every later filter, which is how indels
# stayed invisible even once they could be scored.
VARIANT_ANNOTATION = r'_(?:[A-Z*]\d+_[A-Za-z*]+|\d+(?:-\d+)?_[^;]+):\d+'

LIBRARY_COLUMNS = ["Stripped.Sequence", "Precursor.Charge", "Fragment.Type",
                   "Fragment.Series.Number", "Relative.Intensity", "Decoy"]

# The search cuts after K and R (--cut K*,R*). Only needed to tell whether a truncated protein's
# new C-terminus is novel or is a cleavage site the wild type uses anyway.
CLEAVAGE_RESIDUES = "KR"


def load_variant_sequences(fasta_path):
    """Map full accession -> protein sequence, for wild-type and variant entries alike.

    The FASTA key is the whole accession as it appears in '>sp|P04637_R181_H:9|P53_HUMAN', not
    the base accession, because each variant entry carries its own mutated sequence.
    """
    seqs, key, buf = {}, None, []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if key:
                    seqs[key] = "".join(buf)
                fields = line[1:].split("|")
                key = fields[1] if len(fields) > 1 else line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if key:
        seqs[key] = "".join(buf)
    return seqs


def load_library(parquet_path):
    """Read the DIA-NN library's fragment table, targets only."""
    lib = pq.read_table(parquet_path, columns=LIBRARY_COLUMNS).to_pandas()
    return lib[lib["Decoy"] == 0]


def base_accession(accession):
    """Variant accession -> the wild-type entry it was built from.

    'P35222_23-71_D32_S33del:1' -> 'P35222'.  'Q01196-8_R201_Q:2' -> 'Q01196-8': the isoform
    suffix is part of the accession and must survive, because the wild-type counterpart of an
    isoform-derived variant is that isoform and not the canonical sequence.
    """
    return accession.split("_")[0]


def site_in_peptide(variant_seq, peptide, protein_position):
    """1-based position of the edit within the peptide, or None if the peptide misses the site.

    Substitutions only. Retained for pool_dilution.py; the scoring below works from sequence
    divergence instead and so covers indels too.
    """
    offset = variant_seq.find(peptide)
    if offset < 0:
        return None
    p = protein_position - offset
    return p if 1 <= p <= len(peptide) else None


def _common_prefix(a, b):
    n = min(len(a), len(b))
    i = 0
    while i < n and a[i] == b[i]:
        i += 1
    return i


def _common_suffix(a, b):
    n = min(len(a), len(b))
    i = 0
    while i < n and a[-1 - i] == b[-1 - i]:
        i += 1
    return i


def wt_identical_flanks(variant_seq, wt_seq, peptide):
    """(L, S) -- residues at each end of the peptide whose fragments the wild type also makes.

    None when the peptide cannot testify: absent from the variant sequence, or lying clear of the
    region where variant and wild type diverge, in which case it is a wild-type peptide and no
    fragment of it is evidence of anything.
    """
    offset = variant_seq.find(peptide)
    if offset < 0:
        return None
    n = len(peptide)
    first_diff = _common_prefix(variant_seq, wt_seq)
    last_diff = len(variant_seq) - _common_suffix(variant_seq, wt_seq)

    # A stop-gain truncates the protein, so the variant is a strict prefix of the wild type and
    # there is no divergent residue to span -- the evidence is the new C-terminus itself. The
    # peptide carrying it shares every b-ion with the wild type and no y-ion, since those need a
    # C-terminus the wild type never produces. Any earlier peptide is pure wild-type sequence.
    # This case has to be named: it otherwise looks identical to a peptide sitting upstream of a
    # substitution, which is wild type and must score nothing.
    # Unless the truncation lands straight after a cleavage site: then the wild type cuts there
    # too and makes this very peptide, so not even the C-terminus is novel and nothing is
    # evidence. 27 of the 242 stop-gain entries in the PDX FASTA are of that shape.
    if first_diff == len(variant_seq) < len(wt_seq):
        if offset + n != len(variant_seq) or variant_seq[-1] in CLEAVAGE_RESIDUES:
            return None
        return n, 0

    # Every other class -- substitution, insertion, deletion, delins. The peptide must reach into
    # the divergent region; for a pure deletion that region is empty and the test becomes "the
    # junction falls strictly inside the peptide", which is the same expression.
    if not (offset < last_diff and offset + n > first_diff):
        return None
    L = min(max(first_diff - offset, 0), n)
    S = min(max((offset + n) - last_diff, 0), n)
    return L, S


def spanning_fraction(frags, wt_prefix, wt_suffix):
    """Fraction of a precursor's relative intensity carried by ions the variant changes."""
    spans = [(t == "y" and k > wt_suffix) or (t == "b" and k > wt_prefix)
             for t, k in zip(frags["Fragment.Type"], frags["Fragment.Series.Number"])]
    total = frags["Relative.Intensity"].sum()
    if total <= 0:
        return None
    return float(frags.loc[spans, "Relative.Intensity"].sum() / total)


def peptide_fractions(peptide_list_path, fasta_path, library_path):
    """(stripped sequence, charge) -> spanning fraction, for each variant peptide that resolves.

    Keyed on sequence and charge -- the columns the report and the matrix already carry -- so
    nothing is ever joined on a gene label. That matters: the classifier deliberately relabels a
    shared peptide to whichever paralog the genome confirms while the peptide list keeps the
    accession DIA-NN assigned, so a (gene, protein change) join silently voids every such call.

    Returns (fractions, skipped), skipped being (header, reason) pairs so a peptide that gets no
    value is always accountable rather than quietly absent.
    """
    seqs = load_variant_sequences(fasta_path)
    lib = load_library(library_path)
    fractions, skipped = {}, []

    with open(peptide_list_path) as fh:
        for line in fh:
            header = line.strip()
            if not header:
                continue
            accession, peptide, charge = header.rsplit("_", 2)
            variant_seq = seqs.get(accession)
            if variant_seq is None:
                skipped.append((header, "accession absent from the search FASTA"))
                continue
            wt_seq = seqs.get(base_accession(accession))
            if wt_seq is None:
                skipped.append((header, f"no wild-type entry {base_accession(accession)} "
                                        f"in the search FASTA"))
                continue

            flanks = wt_identical_flanks(variant_seq, wt_seq, peptide)
            if flanks is None:
                skipped.append((header, "peptide does not cover the variant"))
                continue

            charge = int(charge)
            frags = lib[(lib["Stripped.Sequence"] == peptide)
                        & (lib["Precursor.Charge"] == charge)]
            if not len(frags):
                skipped.append((header, "precursor absent from the library"))
                continue
            frac = spanning_fraction(frags, *flanks)
            if frac is None:
                skipped.append((header, "library fragment intensities sum to zero"))
                continue
            # (sequence, charge) is unique in the peptide list; max is a guard, not a policy.
            fractions[(peptide, charge)] = max(frac, fractions.get((peptide, charge), 0.0))

    return fractions, skipped


def failing_sequences(fractions, threshold=FRAGMENT_SPECIFICITY_MIN):
    """Stripped sequences whose every charge state falls below the threshold.

    Aggregated by sequence with max(), matching how the classifier aggregates a call over its
    precursors, so the gated matrix and the classification table cannot disagree about which
    peptides were rejected.
    """
    best = {}
    for (peptide, _charge), frac in fractions.items():
        best[peptide] = max(frac, best.get(peptide, 0.0))
    return {peptide for peptide, frac in best.items() if frac < threshold}
