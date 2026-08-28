# DIA-HMD

Find peptides carrying somatic hotspot mutations in DIA mass-spectrometry data.

DIA-HMD searches your runs with DIA-NN against a proteome that includes hotspot variant
sequences, discards every peptide the normal proteome could have produced just as well, and
reports what is left. Give it the mutations you already know are in the samples and it also
scores each call as a true or false positive.

## Install

Dependencies are Python 3.8+ (pandas, pyarrow) and R 4.0+ (tidyverse, data.table,
RColorBrewer). With conda:

```bash
conda env create -f environment.yml && conda activate diann-pipeline
```

Otherwise `bash install_deps.sh` puts them into an environment you already have. `Dockerfile` and
`apptainer.def` build the same thing as a container; neither contains DIA-NN, so running the
pipeline from one still means mounting DIA-NN in and passing `--runtime native`.

DIA-NN 2.0.2 itself is not bundled. Take the Apptainer image from this repository's releases, or
the native Linux binary from upstream:

```bash
gh release download v1.0 --pattern 'diann-2.0.2.img' --dir .

wget https://github.com/vdemichev/DiaNN/releases/download/2.0/DIA-NN-2.0.2-Academia-Linux.zip
unzip DIA-NN-2.0.2-Academia-Linux.zip
```

Pass whichever you have to `--diann`. The runtime is picked from your PATH: apptainer, then
docker, then the bare binary. Override that with `--runtime apptainer|docker|native`.

## Run it

```bash
bash run.sh --input /path/to/raw_files --diann diann-2.0.2.img --threads 32
```

`--input` is a directory of `.raw.dia` files, all of which go into one search. They can be
symlinks to data held elsewhere; the pipeline resolves them and mounts the real directories into
the container. Results are copied to `results/`.

On a cluster the stages are submitted rather than run, so `run.sh` returns before they finish and
prints the copy command for afterwards.

For an end-to-end test, two COLO205 injections are on Zenodo at
[10.5281/zenodo.19436340](https://doi.org/10.5281/zenodo.19436340) (~12 GB each, too large for
GitHub). [`example/README.md`](example/README.md) has the download and the run, plus a dependency
check that needs no data at all.

## Output

Copied to `results/`:

| File | What it holds |
|---|---|
| `hotspot_peptides.tsv` | intensity matrix for the reported variant peptides |
| `hotspot_peptides_with_canonical.tsv` | the same peptides paired with their wild-type counterparts |
| `hotspot_by_gene.pdf` | mutant vs wild-type scatter plots, one panel per gene |
| `hotspot_by_mutation.pdf` | the same, one panel per mutation |
| `hotspot_detection_classification.tsv` | one row per mutation × sample: q-values, fragment score, runs, TP/FP, filters passed (`--truth` only) |
| `hotspot_detection_classification_summary.txt` | the readable version of that file (`--truth` only) |

One more file is left behind in `scripts/Reports/`: `report_peptidoforms.pr_matrix.strict.tsv`,
DIA-NN's precursor matrix with the rejected variant cells emptied. The plots read it.

## Options

```
  --input DIR         directory of *.raw.dia files (required)
  --diann PATH        DIA-NN image or binary (required)
  --output DIR        where to copy results (default: results/)
  --fasta FILE        search database (default: bundled proteome.fasta)
  --proteome FILE     non-mutated proteome to subtract (default: bundled)
  --runtime RT        apptainer, docker, native, or auto (default: auto)
  --threads N         threads for DIA-NN (default: 4)

  --qvalue-gate Q     per-run Q.Value a call must reach (default: 0.001)
  --fragment-min F    fragment-geometry threshold (default: 0.15)
  --min-replicates N  injections a call must be seen in (default: 1, i.e. off)

  --truth FILE        known mutations, for TP/FP scoring
  --sample-map FILE   which runs are injections of the same sample
  --aliases FILE      truth-table names spelled differently in the run names
  --pools FILE        samples that are mixtures of others
```

To configure by hand instead, edit `scripts/config.sh` and run `bash Complete_pipeline.sh` from
`scripts/`.

## How it works

```
DIA-NN search  →  canonical filter  →  detection gate  →  tables, plots, TP/FP
```

DIA-NN runs twice, a library-free pass to predict a spectral library and a library-guided pass
over your runs. The canonical filter then drops every peptide the non-mutated proteome could also
have produced, treating Ile and Leu as one residue.

### What counts as a detection

What is left is still not evidence of a mutation. The measurement has to be able to tell the
peptide from its normal counterpart, and it has to be a measurement of the sample in question.
Three filters decide that:

- **Q.Value**, per cell, default 0.001. Was this precursor identified in *this* run, at 0.1% run
  FDR?
- **Fragment geometry**, per peptide, default 0.15. Do at least 15% of the library fragment
  intensities come from ions whose mass the mutation changes?
- **Replicates**, per sample, default 1 (off). Was the call seen in at least N injections of the
  sample?

The gate and the classifier read the same thresholds and apply the same rules, so no plot can
show a point the results table rejected.

Replicates are off by default because the pipeline has no way to guess how your injections are
organised: replicates named `S1_1`/`S1_2` and distinct samples named `P_1`/`P_2` look identical
to it. Describe the grouping with `--sample-map`, then raise the value. Asking for
`--min-replicates 2` without a map is an error rather than a run that rejects everything.

[`docs/methods.md`](docs/methods.md) covers where the fragment-geometry rule comes from, how 0.15
was chosen, and what the canonical filter does before any of this.

## Scoring against known mutations

`--truth` turns on the classification stage. It takes a TSV with four columns: `Sample` (matched
against the run-side names), `Gene`, `Protein.Change` (HGVS, e.g. `p.G12V`) and
`Detected.By.DIANN` (`TRUE`/`FALSE`, whether the variant is considered observable at all).

Three further files are for cohorts whose sample names do not line up one-to-one with their run
names. All are optional and all default to the plain case.

| File | Columns | Left out means |
|---|---|---|
| `--sample-map` | `run`, `sample`, optional `run_label` | every run is its own sample |
| `--aliases` | `truth_name`, `run_name` | truth-table names are used as written |
| `--pools` | `pool`, `members` (`*` = any truth mutation, empty = none, else `A;B`) | no sample is a mixture |

Calls in a pool come back as `POOL_TP`/`POOL_FP`, never as `TP`/`FP`: no single genome can
confirm or refute a call made in a mixture.

Both published cohorts ship as worked examples, under `data/truth/` and `data/cohorts/`:

```bash
bash run.sh --input /data/celllines --diann diann-2.0.2.img --threads 32 \
    --truth      data/truth/Table1_hotspotcelllines_stopfree.tsv \
    --sample-map data/cohorts/celllines_sample_map.tsv \
    --aliases    data/cohorts/celllines_aliases.tsv \
    --pools      data/cohorts/celllines_pools.tsv \
    --min-replicates 2
```

`--min-replicates 2` is what the paper used, and it is the only place where this repository's
defaults differ from it. At 1, calls seen in a single injection are admitted: on the same data
that is 22 TP / 4 FP → 24 TP / 7 FP for the PDX cohort, and 6 TP / 1 FP → 6 TP / 2 FP for the
cell lines.

## On a cluster

Stages go through `sbatch` or `qsub` if either is on the PATH, and run one after another in the
current shell if neither is. Nothing needs editing for this.

| Stage | Slots | Memory | Runtime |
|---|---|---|---|
| DIA-NN search | `--threads` | ~6 GB/thread | ~6 h for 20 files at 32 threads |
| Precursor summary | 4 | 32 GB | ~5 min |
| Post-processing | 8 | 32 GB | ~30 min |
| Classification | 4 | 32 GB | ~1 min |
| R analysis | 1 | ~8 GB | ~5 min |

The Python stages ask for several slots to raise the memory ceiling, not to run in parallel:
pyarrow reserves far more virtual memory than it resides, and a 1.5 GB parquet report needs the
whole 32 GB.

SGE starts a held job when the one before it finishes, whatever its exit status, so a failed
search does not stop the chain there. Each stage checks its own inputs first and names the stage
to go and look at. Logs are in `scripts/log/`, one per job.

## What is in `data/`

`fasta/proteome.fasta` is the search database, the reference proteome plus hotspot variant
sequences (23,608 entries). `fasta/proteome_nonmutated.fasta` is its non-mutated subset (20,659),
and is what the canonical filter subtracts. `truth/` and `cohorts/` hold the ground truth, sample
maps, aliases and pool compositions for the cell-line cohort (13 mutations across 6 lines) and
the PDX cohort (64 across 36 samples).

Before you swap the FASTA files out, read the last section of [`docs/methods.md`](docs/methods.md):
the search database deliberately carries no stop-gain variants, and the file being subtracted is
deliberately not a canonical-only proteome.

## Citation

*[to be added on publication]*

## License

MIT. See [LICENSE](LICENSE).
