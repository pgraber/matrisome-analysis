# Environment

Reproducible environment for this compendium. The single source of truth for R package
versions is `../renv.lock`; the `Dockerfile` here wraps it with the OS libraries and the
pipeline runner.

| Component | Version | Pinned in |
|---|---|---|
| R | 4.4.1 | `renv.lock` |
| Bioconductor | 3.19 | `renv.lock` |
| R packages | 148, exact versions (DESeq2, tidyverse, GOplot, circlize, readxl, …) | `renv.lock` |
| Snakemake | 9.9.0 | `Dockerfile` |

## Reproduce locally (fast path)
Local R already matches (4.4.1). From the repo root:
```
R -e "renv::restore(prompt = FALSE)"     # reconstruct the R library from renv.lock
snakemake -s workflow/Snakefile -c1      # run the DAG
```

## Reproduce in the container (clean-room)
```
docker build -t ecm_landscape_phgg -f env/Dockerfile .
docker run --rm -v "$PWD":/project -w /project ecm_landscape_phgg \
  snakemake -s workflow/Snakefile -c1
```

## Notes
- The container base fixes R at 4.4.1. Package versions come from `renv.lock`, so results are
  independent of the base image beyond R's minor version.
- First container build compiles Bioconductor from source (~30–60 min); it is a one-time cost.
- Snakemake is a Python tool; on this machine it lives in the miniforge install, in the
  container it is installed via pip. Either satisfies the pipeline.
