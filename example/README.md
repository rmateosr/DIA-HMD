# Example Data

Example DIA-MS input data is hosted on Zenodo because the files are too large for GitHub (~13 GB).

## Download

```bash
# From the repository root:
wget -O example/24f201_DIA_COLO205.1.raw.dia \
  "https://zenodo.org/records/19436340/files/24f201_DIA_COLO205.1.raw.dia"
wget -O example/24f201_DIA_COLO205.2.raw.dia \
  "https://zenodo.org/records/19436340/files/24f201_DIA_COLO205.2.raw.dia"
```

**Zenodo DOI:** [10.5281/zenodo.19436340](https://doi.org/10.5281/zenodo.19436340)

| File | Size | Description |
|------|------|-------------|
| `24f201_DIA_COLO205.1.raw.dia` | ~12 GB | DIA-MS data from COLO205 cell line (run 1) |
| `24f201_DIA_COLO205.2.raw.dia` | ~12 GB | DIA-MS data from COLO205 cell line (run 2) |

## Running the Example

```bash
# 1. Make sure you have DIA-NN (see the main README for the alternatives)
gh release download v1.0 --pattern 'diann-2.0.2.img' --dir .

# 2. Run the pipeline
bash run.sh --input example/ --diann diann-2.0.2.img --threads 4

# 3. Check results
ls results/
```

Expected outputs in `results/`:

| File | Description |
|------|-------------|
| `hotspot_peptides.tsv` | Intensity matrix for detected hotspot peptides |
| `hotspot_peptides_with_canonical.tsv` | Hotspot peptides with wild-type counterparts |
| `hotspot_by_gene.pdf` | Scatter plots grouped by gene |
| `hotspot_by_mutation.pdf` | Scatter plots grouped by mutation |

`scripts/Reports/report_peptidoforms.pr_matrix.strict.tsv` is also written: DIA-NN's precursor
matrix with the variant cells the detection gate rejected emptied. It is what the two PDFs read.

## Exercising the classification stage

These two files are one cell line injected twice, so they are the case where the replicate
requirement can be switched on — and COLO205's BRAF V600E is in the shipped cell-line truth table,
so the calls can be scored. The shipped sample map names the runs of the full published cohort,
not these two, so write a two-row one for them:

```bash
printf 'run\tsample\n24f201_DIA_COLO205.1\tCOLO205\n24f201_DIA_COLO205.2\tCOLO205\n' \
  > example/colo205_sample_map.tsv

bash run.sh --input example/ --diann diann-2.0.2.img --threads 4 \
  --truth data/truth/Table1_hotspotcelllines_stopfree.tsv \
  --sample-map example/colo205_sample_map.tsv \
  --min-replicates 2
```

The `run` column is the file name with `.raw.dia` removed. That adds
`hotspot_detection_classification.tsv` and `hotspot_detection_classification_summary.txt` to
`results/`; the summary's Section 1 should report BRAF p.V600E as detected in COLO205.

Dropping `--sample-map` and `--min-replicates` also works: each injection is then its own sample,
and the two are reported separately.

## Verifying Your Installation

If you want to check dependencies without downloading the full example data:

```bash
# Check Python deps
python3 -c "import pandas; import pyarrow; print('Python OK')"

# Check the scoring rule and the cohort loaders (need no data at all)
(cd scripts && python3 test_fragment_geometry.py && python3 test_run_samples.py)

# Check R deps
Rscript -e 'library(tidyverse); library(data.table); library(RColorBrewer); cat("R OK\n")'

# Check CLI help
bash run.sh --help

# Check DIA-NN is accessible (no-args prints version banner and exits cleanly)
apptainer exec diann-2.0.2.img /diann-2.0.2/diann-linux
```

## Input Format

The pipeline expects a directory containing one or more `.raw.dia` files from DIA mass
spectrometry experiments. These are typically produced by Thermo Scientific instruments
and converted using appropriate vendor software.
