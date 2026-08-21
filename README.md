# DIANN Paper Pipeline

DIA-NN based proteogenomic pipeline for detecting **somatic mutation (hotspot) peptides** from DIA mass spectrometry data.

## Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/rmateosr/DIA-HMD.git
cd DIA-HMD

# 2. Install dependencies (pick one)
conda env create -f environment.yml && conda activate diann-pipeline
# OR: bash install_deps.sh

# 3. Download external files (not included in the repo due to size)

# DIA-NN 2.0.2 Apptainer image (321 MB) — from the GitHub Release
gh release download v1.0 --pattern 'diann-2.0.2.img' --dir .
# Or download manually from the Releases page:
#   https://github.com/rmateosr/DIA-HMD/releases/tag/v1.0

# Example DIA-MS data (~24 GB total, 2 files) — from Zenodo
wget -O example/24f201_DIA_COLO205.1.raw.dia \
  "https://zenodo.org/records/19436340/files/24f201_DIA_COLO205.1.raw.dia"
wget -O example/24f201_DIA_COLO205.2.raw.dia \
  "https://zenodo.org/records/19436340/files/24f201_DIA_COLO205.2.raw.dia"
# Or download manually from: https://doi.org/10.5281/zenodo.19436340

# 4. Run the pipeline
bash run.sh \
  --input example/ \
  --diann diann-2.0.2.img \
  --threads 4
```

Results will be in `results/`. That run reports every variant peptide that survives the detection
filters. To also have the calls scored against a known truth set, add `--truth` — see
[Ground truth and cohort description](#ground-truth-and-cohort-description).

## What the Pipeline Does

1. **Runs DIA-NN** (two-pass search) against a custom FASTA containing the reference proteome and hotspot variant sequences
2. **Summarises the report** — per-peptide best q-values and the runs each peptide was seen in
3. **Converts** the peptidoform matrix to FASTA format
4. **Filters out canonical peptides** — removes anything a normal protein could also produce, comparing Ile and Leu as one residue because they are isobaric
5. **Applies the detection gate** to the surviving variant peptides — per-run q-value, fragment geometry, and optionally a replicate requirement ([Detection filters](#detection-filters))
6. **Classifies the calls as TP/FP** against a ground truth, if one is supplied
7. **Generates summary tables and plots** — intensity matrices (TSV) and mutant-vs-wild-type scatter plots (PDF)

```
run.sh  (CLI wrapper — parses args, writes config, calls Complete_pipeline.sh)
 └─ Complete_pipeline.sh
     ├─ generate_diann_job.sh          → DIA-NN two-pass search
     ├─ strict_filter.sh               → per-peptide q-value summary (reporting only)
     │    └─ extract_strict_precursors.py
     └─ Post_DIANN_pipeline.sh         → post-processing coordinator
          ├─ awk (inline)              → FASTA conversion
          ├─ filter_canonical_peptides.sh  → canonical peptide removal (I/L collapsed)
          ├─ gate_variant_cells.py         → detection gate → pr_matrix.strict.tsv
          │    └─ fragment_geometry.py     → the fragment-geometry rule
          ├─ classify_job.sh               → TP/FP classification (only with --truth)
          │    └─ classify_hotspot_detections.py
          └─ R hotspot analysis        → hotspot_peptides.tsv + PDFs
```

Both the gate and the classifier apply the same three filters and read the same thresholds, so a
plot cannot show a point the results table rejected.

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
| **Apptainer image (recommended)** | HPC clusters | Download from [this repo's GitHub Release](https://github.com/rmateosr/DIA-HMD/releases/tag/v1.0), or build from the bundled `apptainer.def` |
| **Docker image** | Local machines | `docker pull biocontainers/diann:2.0.2` (check availability) |
| **Native binary** | Any Linux | Download from [DIA-NN releases](https://github.com/vdemichev/DiaNN/releases) |

To download the pre-built Apptainer image from the release:

```bash
# Using GitHub CLI (recommended)
gh release download v1.0 --pattern 'diann-2.0.2.img' --dir .

# Or using curl (GitHub redirects large assets through CDN)
curl -LO "https://github.com/rmateosr/DIA-HMD/releases/download/v1.0/diann-2.0.2.img"

# Or build from the definition file
apptainer build diann-2.0.2.img apptainer.def
```

To use the native binary instead of a container:

```bash
# Download and extract DIA-NN 2.0.2
wget "https://github.com/vdemichev/DiaNN/releases/download/2.0/DIA-NN-2.0.2-Academia-Linux.zip"
unzip DIA-NN-2.0.2-Academia-Linux.zip

# Run pipeline with native binary
bash run.sh --input /path/to/raw_files --diann diann-2.0.2/diann-linux --runtime native
```

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
  --proteome FILE     Non-mutated proteome to subtract (default: bundled)
  --runtime RT        apptainer, docker, native, or auto (default: auto)
  --threads N         Threads for DIA-NN (default: 4)

  --qvalue-gate Q     Per-run Q.Value gate (default: 0.001)
  --fragment-min F    Fragment-geometry threshold (default: 0.15)
  --min-replicates N  Injections a call must be seen in (default: 1 = off)

  --truth FILE        Ground-truth TSV; without it, no TP/FP classification
  --sample-map FILE   Groups injections into samples
  --aliases FILE      Truth-table names spelled differently in the run names
  --pools FILE        Samples that are mixtures of others
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

## Detection filters

A peptide carrying a mutation is only evidence of that mutation if the measurement can tell it
from the normal protein. Three filters decide that, all applied identically by
`gate_variant_cells.py` (which writes the matrix the plots use) and by
`classify_hotspot_detections.py` (which writes the results table).

| Filter | Grain | Default | What it asks |
|---|---|---|---|
| **Q.Value** | per cell | `0.001` | Was this precursor identified in *this* run at 0.1% run FDR? DIA-NN's own output is already at 1%, so this is the gate that decides presence in a sample. `Lib.Q.Value` cannot do it: it carries one value per peptide for the whole study |
| **Fragment geometry** | per peptide | `0.15` | Do at least 15% of the library fragment intensity come from ions whose mass the mutation changes? A fragment that misses the mutated site has the same mass in the mutant as in the wild type, so it is signal the normal protein produces identically and cannot be evidence |
| **Replicate** | per sample | `1` (off) | Was the call seen in at least N injections of the same sample? |

### The fragment-geometry rule

The diagnostic fragment set is derived from where the variant sequence stops matching its
wild-type parent, not from the mutation annotation, so substitutions, insertions, deletions,
delins and stop-gains are all scored by the same arithmetic. With the two proteins sharing a
prefix of `E` residues and a suffix of `F`, and the peptide at offset `o` with length `n`:

```
L = clamp(E - o, 0, n)                        residues still wild type at the N-terminal end
S = clamp((o + n) - (len(variant) - F), 0, n)  ... and at the C-terminal end
b_k differs from wild type iff k > L;  y_k iff k > S
```

The wild-type sequence is taken from the **searched FASTA**, not from a canonical proteome: a
variant built on an isoform compared against the canonical sequence would read every isoform
difference as part of the mutation.

`0.15` is the largest cost-free cut available in both published cohorts — it drops no genomically
confirmed call while taking cell-line precision from 48.0% to 80.0% and PDX from 66.7% to 81.5%.
It was fitted on substitutions; indels are scored against it unchanged. `scripts/test_fragment_geometry.py`
checks the rule on written-out sequences and needs no data:

```bash
cd scripts && python3 test_fragment_geometry.py
```

### The replicate requirement defaults to off

`--min-replicates 1` imposes nothing, because this pipeline cannot know how your injections are
organised — a cohort naming replicates `S1_1`/`S1_2` and one naming distinct samples `P_1`/`P_2`
look the same. To switch it on, describe the grouping with `--sample-map` and raise the value;
`run.sh` refuses `--min-replicates 2` without a map rather than rejecting every call.

**The published cohorts used `--min-replicates 2`.** Relaxing it to 1 admits calls seen in a
single injection: on the same data that is 22 TP / 4 FP → 24 TP / 7 FP for the PDX cohort, and
6 TP / 1 FP → 6 TP / 2 FP for the cell lines. Use 2 with a sample map to reproduce the paper.

## Ground truth and cohort description

`--truth` turns on the TP/FP classification stage. The table needs four columns:

| Column | Meaning |
|---|---|
| `Sample` | Sample name, matched against the run-side names |
| `Gene` | Gene symbol |
| `Protein.Change` | HGVS protein change, e.g. `p.G12V` |
| `Detected.By.DIANN` | `TRUE`/`FALSE` — whether the variant is considered observable at all |

Without `--truth` the stage is skipped and the pipeline still produces the gated matrix, the
hotspot tables and the plots. Three further files describe a cohort whose names do not line up
one-to-one with its runs; all are optional, and absent means the generic case:

| File | Columns | Absent means |
|---|---|---|
| `--sample-map` | `run`, `sample`, optional `run_label` | every run is its own sample |
| `--aliases` | `truth_name`, `run_name` | truth-table names are used as written |
| `--pools` | `pool`, `members` (`*` = any truth mutation, empty = none, else `A;B`) | no sample is a mixture |

Hits in a pool are reported as `POOL_TP`/`POOL_FP` and never as `TP`/`FP`: no single sample's
genome can confirm or refute a call in a mixture. Both published cohorts ship as worked examples
under `data/truth/` and `data/cohorts/`.

## Input Format

### DIA raw files

Directory of `*.raw.dia` DIA mass spectrometry files. All files in the directory are included in the search.

### Example data

Two example DIA-MS files (~12 GB each) are available from Zenodo:

> **DOI:** [10.5281/zenodo.19436340](https://doi.org/10.5281/zenodo.19436340)

```bash
# Download into the example/ directory
wget -O example/24f201_DIA_COLO205.1.raw.dia \
  "https://zenodo.org/records/19436340/files/24f201_DIA_COLO205.1.raw.dia"
wget -O example/24f201_DIA_COLO205.2.raw.dia \
  "https://zenodo.org/records/19436340/files/24f201_DIA_COLO205.2.raw.dia"

# Run the pipeline on them
bash run.sh --input example/ --diann diann-2.0.2.img --threads 4
```

See [`example/README.md`](example/README.md) for details.

## Output

All outputs are written to `results/` (or the directory specified with `--output`):

| File | Description |
|------|-------------|
| `hotspot_peptides.tsv` | Intensity matrix for hotspot (somatic mutation) peptides |
| `hotspot_peptides_with_canonical.tsv` | Hotspot peptides paired with matching wild-type counterparts |
| `hotspot_by_gene.pdf` | Mutant vs wild-type scatter plots, grouped by gene |
| `hotspot_by_mutation.pdf` | Mutant vs wild-type scatter plots, grouped by mutation |
| `hotspot_detection_classification.tsv` | One row per mutation × sample: q-values, fragment score, runs, TP/FP, and which filters it passed. Only with `--truth` |
| `hotspot_detection_classification_summary.txt` | Readable report — recall against the truth table, threshold and filter breakdowns, per-sample detail. Only with `--truth` |

Also written, under `scripts/Reports/` rather than the results directory:

| File | Description |
|------|-------------|
| `report_peptidoforms.pr_matrix.strict.tsv` | DIA-NN's precursor matrix with variant cells the gate rejected emptied. This is what the plots read |
| `strict_precursors_peptide_list.tsv` | Per-peptide summary of best q-values and runs detected |

### Resource usage

| Job | Slots | Memory | Approx. runtime |
|-----|-------|--------|-----------------|
| DIA-NN (two-pass search) | `--threads` value | ~6 GB per thread | ~6 hours (20 DIA files, 32 threads) |
| Precursor summary (Python) | 4 | 32 GB total | ~5 minutes |
| Post-processing (canonical filter + gate) | 8 | 32 GB total | ~30 minutes |
| Classification (Python) | 4 | 32 GB total | ~1 minute |
| R analysis (hotspot) | 1 | ~8 GB | ~5 minutes |

The Python stages take several slots for the memory ceiling, not for parallelism: pyarrow reads a
whole report into memory and reserves far more virtual memory than it resides. A 1.5 GB parquet
report needs the 32 GB; 8 GB was not enough and failed silently.

## Bundled Reference Data

| File | Entries | Description |
|------|---------|-------------|
| `data/fasta/proteome.fasta` | 23,608 | Search database: reference proteome + hotspot variant sequences |
| `data/fasta/proteome_nonmutated.fasta` | 20,659 | The non-mutated subset of the above — canonical sequences plus 15 rescued isoforms — subtracted by the canonical filter |
| `data/truth/Table1_hotspotcelllines_stopfree.tsv` | 13 rows | Ground truth for the published cell-line cohort (6 cell lines) |
| `data/truth/Table_2_HotspotPDX_stopfree.tsv` | 64 rows | Ground truth for the published PDX cohort (36 samples) |
| `data/cohorts/*.tsv` | — | Sample maps, aliases and pool compositions for the two cohorts |

Two things about these files are deliberate and worth knowing before you swap them out.

**The search database carries no stop-gain variants.** A stop-gain entry is a truncated
N-terminal prefix of its wild-type protein, so every peptide it yields is also a canonical peptide
and the canonical filter discards it in every run. Such entries could never produce a reported
call, while still competing in the search, forming `GENE;GENE_var` protein groups, and sitting in
the recall denominator as guaranteed misses. `scripts/make_stopfree_fasta.py` removes them, and
refuses to write if the number it drops is not the number you told it to expect.

**The file the canonical filter subtracts is not a canonical-only proteome.** It is the
non-mutated subset of the database that was actually searched, which includes the 15 isoforms the
library build rescued. Subtracting a canonical-only file instead would report every isoform
wild-type peptide as a variant.

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

## Citation

If you use this pipeline, please cite:

> *[Paper reference to be added upon publication]*

## License

MIT License. See [LICENSE](LICENSE).
