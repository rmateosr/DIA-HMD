# Changelog

## v2.0 — 2026-08-21

The pipeline now decides whether a variant peptide is *evidence* of its mutation, rather than
only whether DIA-NN identified it. Three filters do that, and both the plots and the results
table read the same thresholds, so a figure cannot show a point the results table rejected.

Every change below was checked against the two published cohorts before release. Where a check
is quoted as "identical", it means byte-for-byte against the output the paper was written from.

### Added

- **Fragment-geometry filter** (`scripts/fragment_geometry.py`). A fragment whose mass is the
  same in the mutant as in the wild type cannot be evidence that the variant is present. Which
  fragments those are follows from where the variant sequence stops matching its wild-type
  parent, so substitutions, insertions, deletions, delins and stop-gains are all scored by the
  same arithmetic. A peptide must draw at least 15% of its library fragment intensity from ions
  the mutation changes. On the published cohorts this takes precision from 48.0% to 80.0%
  (cell lines) and 66.7% to 81.5% (PDX) without dropping a single genomically confirmed call.
- **Per-run detection gate** (`scripts/gate_variant_cells.py`). Writes
  `report_peptidoforms.pr_matrix.strict.tsv` with the variant cells that fail the gate emptied.
  The wide matrix holds intensities but no q-values, so this cannot be applied any later.
- **TP/FP classification** (`scripts/classify_hotspot_detections.py`, `scripts/classify_job.sh`),
  enabled by `--truth`. Reports recall against a ground truth, threshold tiers, the filter
  breakdown, and per-mutation and per-sample detail.
- **Optional cohort description** (`scripts/run_samples.py`): `--sample-map`, `--aliases` and
  `--pools`. All optional; absent means the generic case (one sample per run, names as written,
  no mixtures). This replaced two hardcoded cohort-specific copies of the classifier.
- **`scripts/test_fragment_geometry.py`** — checks the rule on written-out sequences, so it runs
  on a fresh clone with no data downloaded.
- **`scripts/make_stopfree_fasta.py`** — builds the search database, and refuses to write if the
  number of stop-gain entries it drops is not the number it was told to expect.
- **Reference data**: `data/truth/` (both cohorts' ground truth) and `data/cohorts/` (sample
  maps, aliases, pool compositions), as worked examples for `--truth`.

### Changed

- **The canonical filter compares Ile and Leu as one residue.** They have identical elemental
  composition and are unresolvable by MS, so a variant peptide that differs from a canonical one
  only by I↔L is not a detection.
- **The canonical filter lays the proteome out one protein per line** instead of flattening it to
  a single string. `grep -F` cannot match past the end of a line, so a peptide can no longer be
  discarded as canonical because it matched a chimera spanning two proteins that no protein
  produces. Records are cut on lines beginning with `>`, not on the character, because twelve
  UniProt descriptions contain one.
- **The bundled search database no longer carries stop-gain variants** (24,050 → 23,608 entries).
  A stop-gain entry is a truncated prefix of its wild-type protein, so every peptide it yields is
  canonical and the filter always discarded it: those entries could only add search-space
  competition and guaranteed misses.
- **`data/fasta/human_canonical_proteome.fasta` → `data/fasta/proteome_nonmutated.fasta`**
  (20,644 → 20,659 entries). The file the canonical filter subtracts must be the non-mutated
  subset of the database that was searched, including the 15 rescued isoforms. A canonical-only
  file reports every isoform wild-type peptide as a variant.
- **`strict_filter.sh` is now reporting only.** Its thresholds equal the 1% FDR DIA-NN already
  applies on output, so it removes nothing; what it writes is a per-peptide summary. The
  filtering that matters moved to the gate, which reads the parquet report and can therefore see
  every individual run a peptide appears in. `Lib.Q.Value` cannot decide presence in a sample: it
  carries one value per peptide for the whole study.
- **`filter_pr_matrix.py` is no longer part of the pipeline** and writes
  `pr_matrix.peptide_filtered.tsv`, so it cannot overwrite the gated matrix.
- **R analysis**: wild-type counterpart genes are resolved through `Protein.Ids` instead of
  DIA-NN's `Genes` field, which names only one of the proteins a peptide is consistent with —
  `IGDFGLATVK` was labelled RAF1 inside the BRAF V600E panel, and unmutated `LVVVGAGGVGK` was
  labelled `KRAS_K117_R:1`. Mutant peptides with no surviving measurement, and wild-type
  counterparts orphaned by that, are dropped with both counts printed.
- **Resources**: the Python stages now take 4 slots for the 32 GB memory ceiling rather than for
  parallelism. pyarrow reserves far more virtual memory than it resides; 8 GB was not enough for
  a 1.5 GB parquet report and the failure was silent.

### Fixed

- R: the second plotting loop indexed the panel's mutation annotation by the outer loop counter,
  reading a different mutation's annotation, or `NA` once the counter passed the subset's row
  count.
- R: `1:length()` over an empty candidate list counted `1, 0`, and a panel with no counterpart
  inherited the previous panel's. Both loops now reset per panel and use `seq_along()`.
- R: `scale_color_manual` aborts rather than recycling when a panel has more series than the
  palette has values, which killed the whole script.

### Notes on reproducing the published cohorts

The replicate requirement defaults to **off** (`--min-replicates 1`), because this pipeline
cannot know how your injections are organised. The published cohorts used `--min-replicates 2`
together with a sample map. Relaxing it to 1 on the same data admits calls seen in a single
injection: 22 TP / 4 FP → 24 TP / 7 FP (PDX) and 6 TP / 1 FP → 6 TP / 2 FP (cell lines).

## v1.0 — 2026-04-08

First public release. See `BUGFIXES_20260406.md` for the portability work that preceded it.
