# ABOUTME: Container image with R + Python deps for the DIANN paper pipeline.
# ABOUTME: Does NOT include DIA-NN — mount your DIA-NN binary or image at runtime.
#
# Build:
#   docker build -t diann-pipeline .
#
# Run (mount your data and DIA-NN image):
#   docker run --rm \
#     -v /path/to/raw/files:/data/input \
#     -v /path/to/output:/data/output \
#     -v /path/to/diann-linux:/opt/diann/diann-linux \
#     diann-pipeline \
#     --input /data/input --output /data/output \
#     --diann /opt/diann/diann-linux --runtime native

FROM rocker/r-ver:4.4.3

LABEL maintainer="Raul N. Mateos"
LABEL description="DIANN paper pipeline — proteogenomic hotspot peptide detection"

# System deps for R packages (tidyverse needs libcurl, libxml2, etc.)
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Python deps
RUN pip3 install --no-cache-dir --break-system-packages \
    pandas==2.2.* \
    pyarrow==17.*

# R deps
RUN Rscript -e 'install.packages(c("tidyverse", "RColorBrewer", "data.table"), repos="https://cloud.r-project.org", quiet=TRUE)'

# Copy pipeline
COPY . /opt/pipeline
WORKDIR /opt/pipeline
RUN chmod +x run.sh scripts/*.sh

ENV PATH="/opt/pipeline:${PATH}"

ENTRYPOINT ["bash", "run.sh"]
