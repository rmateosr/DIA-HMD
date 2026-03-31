# DIANN Paper Pipeline

DIA-NN based proteogenomic pipeline for detecting **somatic mutation (hotspot) peptides** from DIA mass spectrometry data.

## Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/rmateosr/DIANN_pipeline.git
cd DIANN_pipeline

# 2. Install dependencies (pick one)
conda env create -f environment.yml && conda activate diann-pipeline
# OR: bash install_deps.sh

# 3. Run the pipeline
bash run.sh \
  --input /path/to/your/raw_dia_files \
  --diann /path/to/diann-2.0.2.img \
  --threads 4
```

Results will be in `results/`.

## What the Pipeline Does

1. **Runs DIA-NN** (two-pass search) against a custom FASTA containing the reference proteome and hotspot variant sequences
2. **Applies strict q-value filtering** — filters precursors by Lib.Q.Value (0.1% FDR) for statistically correct MBR results
3. **Converts** the filtered peptidoform matrix to FASTA format
4. **Filters out canonical peptides** — removes exact matches to the UniProt human reference proteome
5. **Identifies non-canonical hotspot hits** among the remaining peptides
6. **Generates summary tables and plots** — intensity matrices (TSV) and mutant-vs-wild-type scatter plots (PDF)

```
run.sh  (CLI wrapper — parses args, writes config, calls Complete_pipeline.sh)
 └─ Complete_pipeline.sh
     ├─ generate_diann_job.sh          → DIA-NN two-pass search
     ├─ strict_filter.sh               → q-value filtering (Python)
     │    ├─ extract_strict_precursors.py  → Lib.Q.Value filter on parquet
     │    └─ filter_pr_matrix.py           → subset pr_matrix to strict peptides
     └─ Post_DIANN_pipeline.sh         → post-processing coordinator
          ├─ awk (inline)              → FASTA conversion
          ├─ filter_canonical_peptides.sh  → canonical peptide removal
          └─ R hotspot analysis        → hotspot_peptides.tsv + PDFs
```

## Installation

### Option A: Conda (recommended)

The simplest way to install all dependencies. Works on Linux and macOS.

```bash
conda env create -f environment.yml
conda activate diann-pipeline
```

This installs Python 3.12, R 4.4, and all required packages.

### Option B: Manual installation

If you don't use conda:

```bash
bash install_deps.sh
```

This installs via `pip` (Python packages) and `install.packages()` (R packages). Requires Python >= 3.8 and R >= 4.0 already on your PATH.

### Option C: Docker container

For a fully reproducible environment:

```bash
docker build -t diann-pipeline .
```

See [Running with Docker](#running-with-docker) below.

### DIA-NN installation

This pipeline requires [DIA-NN 2.0.2](https://github.com/vdemichev/DiaNN). You need one of:

| Method | Best for | How to get |
|--------|----------|------------|
| **Apptainer image** | HPC clusters | Build from DIA-NN releases or pull from a registry |
| **Docker image** | Local machines | `docker pull biocontainers/diann:2.0.2` (check availability) |
| **Native binary** | Any Linux | Download from [DIA-NN releases](https://github.com/vdemichev/DiaNN/releases) |

## Usage

### Basic usage (CLI wrapper)

```bash
bash run.sh --input /path/to/raw_files --diann /path/to/diann-2.0.2.img
```

The `run.sh` wrapper handles configuration, runs the pipeline, and copies results to `results/`.

```
Options:
  --input DIR         Directory containing *.raw.dia files (required)
  --diann PATH        DIA-NN image or binary path (required)
  --output DIR        Where to copy results (default: results/)
  --fasta FILE        Custom FASTA (default: bundled proteome.fasta)
  --proteome FILE     Canonical proteome (default: bundled)
  --runtime RT        apptainer, docker, native, or auto (default: auto)
  --threads N         Threads for DIA-NN (default: 4)
  --help              Show help
```

### Container runtime auto-detection

The pipeline auto-detects how to run DIA-NN:

| Runtime | Detection | `--diann` value |
|---------|-----------|-----------------|
| **Apptainer** | `apptainer` on PATH | Path to `.img` or `.sif` file |
| **Docker** | `docker` on PATH | Image name (e.g., `biocontainers/diann:2.0.2`) |
| **Native** | fallback | Path to `diann-linux` binary |

Override with `--runtime apptainer|docker|native`.

### Running with Docker

If you built the Docker image (Option C), the R/Python environment is inside the container. You still need DIA-NN separately:

```bash
# With a native DIA-NN binary mounted into the container
docker run --rm \
  -v /path/to/raw_files:/data/input \
  -v /path/to/output:/data/output \
  -v /path/to/diann-linux:/opt/diann/diann-linux \
  diann-pipeline \
  --input /data/input --output /data/output \
  --diann /opt/diann/diann-linux --runtime native
```

### Advanced: direct script execution

If you prefer to configure manually instead of using `run.sh`:

```bash
# 1. Edit scripts/config.sh (set SAMPLE_DIR, DIANN_IMG, etc.)
vi scripts/config.sh

# 2. Run from the scripts directory
cd scripts/
bash Complete_pipeline.sh
```

### Scheduler support

The pipeline automatically detects and uses whichever scheduler is available:

| Scheduler | Detection | How it works |
|-----------|-----------|--------------|
| **SLURM** | `sbatch` on PATH | Jobs submitted via `sbatch` with `--dependency=afterok` chains |
| **SGE** | `qsub` on PATH | Jobs submitted via `qsub` with `-hold_jid` chains |
| **None** | fallback | Jobs run sequentially in the current shell (suitable for local machines) |

No manual editing of scheduler directives is needed. On a local machine without a scheduler, all steps run sequentially.

## Input Format

### DIA raw files

Directory of `*.raw.dia` DIA mass spectrometry files. All files in the directory are included in the search.

## Output

All outputs are written to `results/` (or the directory specified with `--output`):

| File | Description |
|------|-------------|
| `hotspot_peptides.tsv` | Intensity matrix for hotspot (somatic mutation) peptides |
| `hotspot_peptides_with_canonical.tsv` | Hotspot peptides paired with matching wild-type counterparts |
| `hotspot_by_gene.pdf` | Mutant vs wild-type scatter plots, grouped by gene |
| `hotspot_by_mutation.pdf` | Mutant vs wild-type scatter plots, grouped by mutation |

### Resource usage

| Job | Threads | Memory | Approx. runtime |
|-----|---------|--------|-----------------|
| DIA-NN (two-pass search) | `--threads` value | ~4 GB per thread | ~6 hours (20 DIA files, 32 threads) |
| Strict filtering (Python) | 1 | ~8 GB | ~5 minutes |
| Post-processing (filtering + R) | 8 | ~32 GB | ~30 minutes |
| R analysis (hotspot) | 1 | ~8 GB | ~5 minutes |

## Bundled Reference Data

| File | Size | Description |
|------|------|-------------|
| `data/fasta/proteome.fasta` | 16 MB | Combined reference + hotspot variant sequences |
| `data/fasta/human_canonical_proteome.fasta` | 13 MB | UniProt canonical human proteome (for filtering) |

## Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| [DIA-NN](https://github.com/vdemichev/DiaNN) | 2.0.2 | DIA proteomics search engine |
| Python | >= 3.8 | Strict q-value filtering |
| R | >= 4.0 | Hotspot analysis and visualization |
| pandas | >= 2.0 | Tabular data processing (Python) |
| pyarrow | >= 12.0 | Parquet file reading (Python) |
| tidyverse | >= 2.0 | Data wrangling and plotting (R) |
| RColorBrewer | >= 1.1 | Color palettes (R) |
| data.table | >= 1.14 | Fast file reading (R) |
| Apptainer or Docker | any | Container runtime for DIA-NN (optional if using native binary) |

## Gene Fusion Analysis

A version of this pipeline that includes gene-fusion peptide detection is preserved in the [`with-fusions`](https://github.com/rmateosr/DIANN_pipeline/tree/with-fusions) branch (tag: `v1.0-with-fusions`).

## Citation

If you use this pipeline, please cite:

> *[Paper reference to be added upon publication]*

## License

MIT License. See [LICENSE](LICENSE).
