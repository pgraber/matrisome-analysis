# Environment for ecm_landscape_phgg (thesis Chapter 3 — ECM/matrisome transcriptional landscape of pHGG)
# Captures the exact analysis environment: R 4.4.1 + Bioconductor 3.19 + 148 pinned
# packages (renv.lock) + Snakemake 9.9.0. Package versions are pinned by renv.lock, so
# reproducibility does not depend on this base image beyond the R minor version.
#
# Build (from repo root):   docker build -t ecm_landscape_phgg -f env/Dockerfile .
# Run    (from repo root):  docker run --rm -v "$PWD":/project -w /project ecm_landscape_phgg \
#                             snakemake -s workflow/Snakefile -c1
#
# NOTE: first build compiles Bioconductor packages from source (~30–60 min). This is a
# one-time cost for an exact-version restore; subsequent builds reuse the renv cache layer.

FROM rocker/r-ver:4.4.1

# --- system libraries needed by tidyverse, DESeq2, GOplot/circlize, readxl, fonts ---
RUN apt-get update && apt-get install -y --no-install-recommends \
      libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
      libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
      libharfbuzz-dev libfribidi-dev libcairo2-dev \
      zlib1g-dev libbz2-dev liblzma-dev \
      python3 python3-pip \
    && rm -rf /var/lib/apt/lists/*

# --- Snakemake (pipeline runner) — version-pinned to match the authoring environment ---
RUN pip3 install --no-cache-dir --break-system-packages snakemake==9.9.0

# --- R packages: restore EXACT versions from the lockfile ---
WORKDIR /project
RUN R -e "install.packages('renv', repos='https://packagemanager.posit.co/cran/__linux__/jammy/latest')"
COPY renv.lock renv.lock
RUN R -e "renv::restore(lockfile = 'renv.lock', prompt = FALSE)"

CMD ["bash"]
