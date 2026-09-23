# Matrisome Analysis

A Snakemake workflow for differential expression and extracellular matrix (matrisome) analysis of
bulk RNA-seq data, comparing paediatric tumour samples against normal brain reference samples.

## What it does

The workflow runs end to end from raw count matrices to publication figures:

1. **Differential expression** (`src/differential_expression.R`). DESeq2 on tumour versus normal
   brain counts, producing the full and significance-filtered result tables.
2. **Functional enrichment** (`src/enrichment_analysis.R`). g:Profiler results rendered as GOplot
   chord diagrams for molecular function, biological process and cellular component.
3. **Matrisome annotation** (`src/annotate_matrisome.R`). Annotates differentially expressed genes
   against the matrisome classification, separating core matrisome from matrisome-associated genes.
4. **Figures** (`src/deg_barplot.R`, `src/matrisome_dotplot.R`, `src/core_matrisome_dotplot.R`,
   `src/pvalue_dotplot.R`). Barplots of core matrisome differential expression and dotplots of
   effect size and significance across matrisome categories.

Every step is declared as a Snakemake rule with explicit inputs and outputs, so the dependency graph
is the pipeline and any step can be rerun in isolation.

## Running it

```
snakemake --cores 4
```

`rule all` declares the full target set. Outputs are written to `output/`.

## Input data

Input count matrices and metadata are not included in this repository. The workflow expects:

- tumour RNA-seq counts, from a controlled-access cohort
- normal brain RNA-seq counts, from a public reference resource
- a sample metadata table
- g:Profiler enrichment exports

Paths are declared at the top of each rule in the `Snakefile`. Controlled-access data is never
committed here.

## Author

Philipp Graber. [ORCID 0000-0002-3157-1434](https://orcid.org/0000-0002-3157-1434)
