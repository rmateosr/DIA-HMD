# Example Data

Example input data is not included because DIA mass spectrometry files use a proprietary
binary format (`.raw.dia`) that cannot be synthesized.

## Verifying Your Installation

After installing dependencies (conda or manual), verify everything works:

```bash
# Check Python deps
python3 -c "import pandas; import pyarrow; print('Python OK')"

# Check R deps
Rscript -e 'library(tidyverse); library(data.table); library(RColorBrewer); cat("R OK\n")'

# Check CLI help
bash run.sh --help

# Check DIA-NN is accessible (replace with your path)
apptainer exec /path/to/diann-2.0.2.img /diann-2.0.2/diann-linux --help
# or for Docker:
docker run --rm biocontainers/diann:2.0.2 /diann-2.0.2/diann-linux --help
# or for native binary:
/path/to/diann-linux --help
```

## Input Format

The pipeline expects a directory containing one or more `.raw.dia` files from DIA mass
spectrometry experiments. These are typically produced by Thermo Scientific instruments
and converted using appropriate vendor software.
