# How a detection is decided

Everything here is about one question: when does an identified variant peptide count as evidence
that the mutation is present in that sample? DIA-NN answers a different question — whether the
peptide was identified at 1% FDR across the study — and the gap between the two is where the
false positives live. This page describes what closes it.

## The canonical filter

The search database contains variant sequences and normal ones, so most of what comes back is
ordinary proteome. `filter_canonical_peptides.sh` removes every peptide that the non-mutated
subset of the searched database could also have produced.

Two details matter more than they look:

**Isoleucine and leucine are compared as one residue.** They have identical elemental
composition and mass spectrometry cannot separate them, so a variant peptide differing from a
canonical one only by I↔L is not a detection. Both operands get `tr 'I' 'L'` before the
comparison.

**The proteome is laid out one protein per line, not flattened into a single string.** `grep -F`
will not match past the end of a line, so a peptide is only ever compared against a sequence
that actually exists. Flattening lets a peptide straddling two proteins match a chimera nothing
produces, and be thrown away as canonical on that evidence. Records are cut on lines beginning
with `>` rather than on the character itself, because twelve UniProt descriptions contain one
(`DNA dC->dU-editing enzyme`) and splitting there breaks those proteins in two.

## The three filters

What survives the canonical filter goes through three more checks. `gate_variant_cells.py`
applies them to the precursor matrix the plots read, `classify_hotspot_detections.py` applies
them to the results table, and both import their thresholds from the same place.

**Q.Value, per cell, default 0.001.** Was this precursor identified in *this* run? DIA-NN's
output is already filtered at 1%, so the study-wide numbers cannot distinguish a sample where
the peptide is present from one where it is not. `Lib.Q.Value` cannot do it either: it carries
one value per peptide for the whole study. Only `Q.Value`, recomputed inside each run and
available in the parquet report, is per-run — which is also why the gate has to be applied
before the wide matrix, which holds intensities but no q-values.

**Fragment geometry, per peptide, default 0.15.** Derived below.

**Replicates, per sample, default 1 (off).** Was the call seen in at least N injections of the
same sample? Off by default because the pipeline cannot know your replicate structure; the
published cohorts used 2, with a sample map. The gate applies this per precursor while the
classifier counts runs per call across all of that call's precursors. The two agree whenever
some single precursor already reaches the required number of runs, which held for every reported
call in both cohorts, and where they could differ the script says so rather than diverging
quietly.

## The fragment-geometry rule

A fragment whose mass is the same in the mutant as in the wild type cannot be evidence that the
variant is present: it is signal the normal protein produces identically. Which fragments those
are follows from where the variant sequence stops matching its wild-type parent. If the two
proteins share a prefix of `E` residues and a suffix of `F`, and the peptide sits at offset `o`
with length `n`:

```
L = clamp(E - o, 0, n)                          residues still wild type at the N-terminal end
S = clamp((o + n) - (len(variant) - F), 0, n)   and at the C-terminal end

b_k differs from wild type iff k > L
y_k differs from wild type iff k > S
```

Ion type and series number are all this needs — no residue masses, no modifications, no
isotopes. A peptide passes if the fragments that differ carry at least `--fragment-min` of the
library fragment intensity.

Taking `L` and `S` from the sequences rather than from the accession's mutation annotation is
what lets substitutions, insertions, deletions, delins and stop-gains be scored by the same
arithmetic. For a substitution at peptide position `p` it collapses to `L = p-1`, `S = n-p`,
the older `b_k >= p` / `y_k >= n-p+1` form, and reproduces every score that form produced.

The wild-type sequence has to come from the searched FASTA and not from a canonical proteome. A
variant built on an isoform (`Q01196-8`) compared against the canonical sequence reads every
isoform difference as part of the mutation; that alone moved `ITVDGPQEPR` from 0.967 to 0.985.

Two limits worth knowing. The score is a property of the library entry, so there is one value
per precursor for the whole study — it can remove peptides incapable of proving their variant,
but it can never say why a peptide is right in sample A and wrong in sample B. And charge states
are separate library rows, counted separately, so the denominator is a sum of ratios rather than
a physical intensity; y7 at 1+ and 2+ can together account for ~46% of one peptide's total.
Collapsing charge states into one entry per (type, series) would change every value and has not
been tested.

A mass-arithmetic version of the rule was written and rejected. It needs a residue-mass table
and fixed-modification handling, and changes no result: its only extra catch is I↔L, which the
canonical filter has already removed.

### Where 0.15 comes from

0.15 is the largest threshold that costs nothing on the two published cohorts. It drops no
genomically confirmed call while taking precision from 48.0% to 80.0% on the cell lines and from
66.7% to 81.5% on the PDX samples. What caps it is the class of mutations that removes a
cleavage site (reference residue K or R): the mutant peptide fuses two wild-type peptides, so
only the weak long b-ions are diagnostic, and the lowest true call of that class, PIK3R1 p.K567E,
sits at 17.9%.

It was chosen on the same cohorts those precision figures describe, and on substitutions only.
Indels are scored against it unchanged. Treat the numbers as a description of these two datasets
rather than as a held-out estimate of what the threshold will do on yours.

`scripts/test_fragment_geometry.py` checks the rule against written-out sequences and
`scripts/test_run_samples.py` checks the cohort loaders. Neither needs any data:

```bash
cd scripts && python3 test_fragment_geometry.py && python3 test_run_samples.py
```

## The two bundled FASTA files

**The search database carries no stop-gain variants.** A stop-gain entry is a truncated
N-terminal prefix of its wild-type protein, so every peptide it can yield is also a canonical
peptide and the canonical filter discards it in every run. Such entries could never produce a
reported call, while still competing in the search, forming `GENE;GENE_var` protein groups and
sitting in the recall denominator as guaranteed misses. `scripts/make_stopfree_fasta.py` removes
them, and refuses to write if the number it drops is not the number you told it to expect.

**The file the canonical filter subtracts is not a canonical-only proteome.** It is the
non-mutated subset of the database that was actually searched, which includes the 15 isoforms
the library build rescued. Subtract a canonical-only file instead and every wild-type peptide
from those isoforms gets reported as a variant.
