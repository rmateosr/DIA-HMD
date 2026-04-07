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
# 1. Make sure you have the DIA-NN image (download from the GitHub Release or build it)
gh release download v1.0 --pattern 'diann-2.0.2.img' --dir .
# Or: apptainer build diann-2.0.2.img apptainer.def

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

## Verifying Your Installation

If you want to check dependencies without downloading the full example data:

```bash
# Check Python deps
python3 -c "import pandas; import pyarrow; print('Python OK')"

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
