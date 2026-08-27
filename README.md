# DIA-HMD

Finding peptides that carry somatic hotspot mutations in DIA mass spectrometry data.

DIA-HMD runs a DIA-NN search against a proteome that includes hotspot variant sequences,
discards every peptide the normal proteome could have produced just as well, and reports what
is left. Hand it the mutations you already know are in your samples and it will also score each
call as a true or false positive.

## Installing

Python and R dependencies, either with conda

```bash
conda env create -f environment.yml
conda activate diann-pipeline
```

or into an environment you already have (Python >= 3.8, R >= 4.0)

```bash
bash install_deps.sh
```

Either route installs pandas and pyarrow on the Python side, tidyverse, data.table and
RColorBrewer on the R side. `Dockerfile` and `apptainer.def` build the same environment as a
container if that suits your site better. Neither of them contains DIA-NN, so if you run the
pipeline from the container, DIA-NN still has to be mounted in from outside:

```bash
docker run --rm \
  -v /path/to/raw_files:/data/input \
  -v /path/to/output:/data/output \
  -v /path/to/diann-linux:/opt/diann/diann-linux \
  diann-pipeline \
  --input /data/input --output /data/output \
  --diann /opt/diann/diann-linux --runtime native
```

## Getting DIA-NN

The pipeline needs DIA-NN 2.0.2, which is not bundled. Two ways to get it:

```bash
# an Apptainer image, from this repository's release
gh release download v1.0 --pattern 'diann-2.0.2.img' --dir .

# or the native Linux binary, from upstream
wget https://github.com/vdemichev/DiaNN/releases/download/2.0/DIA-NN-2.0.2-Academia-Linux.zip
unzip DIA-NN-2.0.2-Academia-Linux.zip
```

Pass whichever you have to `--diann`. How it gets executed follows from what is on your PATH:
apptainer, then docker, then the bare binary. Use `--runtime apptainer|docker|native` to
decide yourself.

## Running it

```bash
bash run.sh --input /path/to/raw_files --diann diann-2.0.2.img --threads 32
```

`--input` is a directory of `.raw.dia` files, all of which go into one search. They can be
symlinks to data held elsewhere; the pipeline resolves them and mounts the real directories
into the container.

Results are copied to `results/`. On a cluster the stages are submitted rather than run, so
`run.sh` returns before they finish and prints the copy command for afterwards.

### On two example files

Two COLO205 injections are on Zenodo ([10.5281/zenodo.19436340](https://doi.org/10.5281/zenodo.19436340),
~12 GB each), too large for GitHub:

```bash
wget -O example/24f201_DIA_COLO205.1.raw.dia \
  "https://zenodo.org/records/19436340/files/24f201_DIA_COLO205.1.raw.dia"
wget -O example/24f201_DIA_COLO205.2.raw.dia \
  "https://zenodo.org/records/19436340/files/24f201_DIA_COLO205.2.raw.dia"

bash run.sh --input example/ --diann diann-2.0.2.img --threads 4
```

COLO205's BRAF V600E is in the shipped cell-line truth table, so this run can also be scored —
see [`example/README.md`](example/README.md) for that and for a dependency check that needs no
data at all.

## What comes out

| File | What it holds |
|---|---|
| `hotspot_peptides.tsv` | intensity matrix for the reported variant peptides |
| `hotspot_peptides_with_canonical.tsv` | the same peptides paired with their wild-type counterparts |
| `hotspot_by_gene.pdf` | mutant vs wild-type scatter plots, one panel per gene |
| `hotspot_by_mutation.pdf` | the same, one panel per mutation |
| `hotspot_detection_classification.tsv` | one row per mutation × sample: q-values, fragment score, runs, TP/FP, filters passed (needs `--truth`) |
| `hotspot_detection_classification_summary.txt` | the readable version of that file (needs `--truth`) |

Two more are left in `scripts/Reports/` rather than copied: `report_peptidoforms.pr_matrix.strict.tsv`,
DIA-NN's precursor matrix with the rejected variant cells emptied, which is what the plots read,
and `strict_precursors_peptide_list.tsv`, a per-peptide summary of best q-values and the runs each
peptide was seen in.

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
the `scripts/` directory.

## How it runs

`run.sh` checks the arguments, writes `scripts/config.run.sh` and hands over to
`Complete_pipeline.sh`, which chains three jobs:

1. `generate_diann_job.sh` writes the search job. DIA-NN runs twice: a library-free pass over the
   search database to predict a spectral library, then a library-guided pass over your runs.
2. `strict_filter.sh` summarises the report as best q-values and runs per peptide. This is
   reporting only, and filters nothing.
3. `Post_DIANN_pipeline.sh` does the rest. It converts the peptidoform matrix to FASTA, drops
   canonical peptides with `filter_canonical_peptides.sh`, applies the detection gate with
   `gate_variant_cells.py`, and submits the TP/FP classification (`classify_job.sh`, only with
   `--truth`) and the R analysis that writes the tables and plots.

## What counts as a detection

An identified variant peptide is not by itself evidence of the mutation. The measurement has to
be able to tell the peptide from its normal counterpart, and it has to be a measurement of the
sample in question. Three filters decide that:

| Filter | Grain | Default | Question |
|---|---|---|---|
| Q.Value | per cell | 0.001 | was this precursor identified in *this* run, at 0.1% run FDR? |
| Fragment geometry | per peptide | 0.15 | do at least 15% of the library fragment intensities come from ions whose mass the mutation changes? |
| Replicates | per sample | 1 (off) | was the call seen in at least N injections of the sample? |

The gate that writes the matrix and the classifier that writes the results table read the same
thresholds and apply the same rules, so no plot can show a point the table rejected.
[`docs/methods.md`](docs/methods.md) explains where the fragment-geometry rule comes from, how
0.15 was chosen and what the canonical filter does before any of this.

The replicate requirement is off by default because the pipeline has no way to guess how your
injections are organised: a cohort naming replicates `S1_1`/`S1_2` and one naming distinct
samples `P_1`/`P_2` look identical to it. Describe the grouping with `--sample-map` and then
raise the value. Asking for `--min-replicates 2` without a map is an error rather than a run
that rejects everything.

## Scoring against known mutations

`--truth` takes a TSV with four columns and turns on the classification stage:

| Column | Meaning |
|---|---|
| `Sample` | sample name, matched against the run-side names |
| `Gene` | gene symbol |
| `Protein.Change` | HGVS protein change, e.g. `p.G12V` |
| `Detected.By.DIANN` | `TRUE`/`FALSE`, whether the variant is considered observable at all |

Three further files are for cohorts whose sample names do not line up one-to-one with their run
names. All are optional and all default to the plain case.

| File | Columns | Left out means |
|---|---|---|
| `--sample-map` | `run`, `sample`, optional `run_label` | every run is its own sample |
| `--aliases` | `truth_name`, `run_name` | truth-table names are used as written |
| `--pools` | `pool`, `members` (`*` = any truth mutation, empty = none, else `A;B`) | no sample is a mixture |

Calls in a pool come back as `POOL_TP`/`POOL_FP`, never as `TP`/`FP`: no single genome can
confirm or refute a call made in a mixture.

## Reproducing the published cohorts

Both cohorts ship as worked examples, under `data/truth/` and `data/cohorts/`:

```bash
bash run.sh --input /data/celllines --diann diann-2.0.2.img --threads 32 \
    --truth      data/truth/Table1_hotspotcelllines_stopfree.tsv \
    --sample-map data/cohorts/celllines_sample_map.tsv \
    --aliases    data/cohorts/celllines_aliases.tsv \
    --pools      data/cohorts/celllines_pools.tsv \
    --min-replicates 2
```

`--min-replicates 2` is what the paper used, and it is the only place where this repository's
defaults differ from it. At 1, calls seen in a single injection are admitted: on the same data that is
22 TP / 4 FP → 24 TP / 7 FP for the PDX cohort and 6 TP / 1 FP → 6 TP / 2 FP for the cell lines.

## On a cluster

Stages are submitted through `sbatch` or `qsub` if either is on the PATH, and run one after
another in the current shell if neither is. Nothing needs editing for this.

| Stage | Slots | Memory | Runtime |
|---|---|---|---|
| DIA-NN search | `--threads` | ~6 GB/thread | ~6 h for 20 files at 32 threads |
| Precursor summary | 4 | 32 GB | ~5 min |
| Post-processing | 8 | 32 GB | ~30 min |
| Classification | 4 | 32 GB | ~1 min |
| R analysis | 1 | ~8 GB | ~5 min |

The Python stages ask for several slots to raise the memory ceiling, not to run in parallel:
pyarrow reads a whole report into memory and reserves far more virtual memory than it resides.
A 1.5 GB parquet report needs the 32 GB. With 8 GB it died silently.

SGE's `-hold_jid` releases a job when the one before it finishes, whatever its exit status,
unlike SLURM's `--dependency=afterok`. A failed search therefore does not stop the chain on SGE,
so every stage checks its own inputs first and names the stage to go and look at:

```
ERROR: required input is missing or empty: Reports/report_peptidoforms.parquet
ERROR: the DIANN stage did not produce what this stage needs.
       Look at log/DIANN.* for why, then rerun from there.
```

Logs are in `scripts/log/`, one per job.

## What is in `data/`

| File | Entries | What it is |
|---|---|---|
| `fasta/proteome.fasta` | 23,608 | the search database: reference proteome plus hotspot variant sequences |
| `fasta/proteome_nonmutated.fasta` | 20,659 | its non-mutated subset, which is what the canonical filter subtracts |
| `truth/Table1_hotspotcelllines_stopfree.tsv` | 13 | known mutations for the cell-line cohort (6 cell lines) |
| `truth/Table_2_HotspotPDX_stopfree.tsv` | 64 | known mutations for the PDX cohort (36 samples) |
| `cohorts/*.tsv` | 4 | sample maps, aliases and pool compositions for those two cohorts |

Before you swap the FASTA files out, read the last section of [`docs/methods.md`](docs/methods.md):
the search database deliberately carries no stop-gain variants, and the file being subtracted is
deliberately not a canonical-only proteome.

## Citation

*[to be added on publication]*

## License

MIT. See [LICENSE](LICENSE).
