#!/bin/bash
# ABOUTME: Installs R and Python dependencies for the DIANN paper pipeline.
# ABOUTME: Run once before using the pipeline. Requires R >= 4.x and Python >= 3.8.
set -euo pipefail

echo "=== DIANN Paper Pipeline — Dependency Installer ==="
echo ""

# --- Python packages ---
echo "Installing Python packages..."
pip install pandas pyarrow
echo "Python packages installed."
echo ""

# --- R packages ---
echo "Installing R packages..."
Rscript -e '
cran_pkgs <- c("tidyverse", "RColorBrewer", "data.table")

# Install CRAN packages
missing_cran <- cran_pkgs[!cran_pkgs %in% installed.packages()[, "Package"]]
if (length(missing_cran) > 0) {
  cat("Installing CRAN packages:", paste(missing_cran, collapse=", "), "\n")
  install.packages(missing_cran, repos="https://cloud.r-project.org", quiet=TRUE)
} else {
  cat("All CRAN packages already installed.\n")
}

cat("\nAll R dependencies installed successfully.\n")
'

echo ""
echo "=== Done. All dependencies installed. ==="
