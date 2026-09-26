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
conda env create -f environment.yml && conda activate dia-hmd
```

Otherwise `bash install_deps.sh` puts them into an environment you already have. `Dockerfile` and
`apptainer.def` build the same thing as a container; neither contains DIA-NN, so running the
pipeline from one still means mounting DIA-NN in and passing `--runtime native`.

DIA-NN 2.0.2 itself is not bundled. Take the Apptainer image from this repository's releases, or
the native Linux binary from upstream:

```bash
wget https://github.com/rmateosr/DIA-HMD/releases/download/v1.0/diann-2.0.2.img

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

## Scoring against known mutations

`--truth` takes a TSV with `Sample`, `Gene`, `Protein.Change` and `Detected.By.DIANN` columns —
one row per mutation you already know is in a sample. With it, every call is reported as a true
or false positive; without it that stage is skipped. `--sample-map` (`run`, `sample`) groups
injections of the same sample, which is what `--min-replicates` counts over.
`example/colo205_truth.tsv` and `example/colo205_sample_map.tsv` are worked ones.

## Citation

*[to be added on publication]*

## License

MIT. See [LICENSE](LICENSE).
