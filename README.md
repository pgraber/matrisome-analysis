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

## Validation

Every step asserts the shape of what it produced, immediately after producing it, against a value
known from outside the code that produced it. Checks print to the terminal as the workflow runs and
append to `output/checks.log` and `output/checks.tsv`. A failed check stops the run, because nothing
downstream should execute on data of the wrong shape.

What is asserted: the tumour sample count against the metadata sheet, gene count after the
tumour-to-normal join against the identifier intersection computed before the join, uniqueness of
gene IDs, non-negativity of counts, p-values within their defined range, and that matrisome
annotation neither drops nor duplicates a gene. The helper is `scripts/validate.R`.

`output/checks.log` and `output/checks.tsv` are the only outputs committed. Everything else in
`output/` is gitignored.

## Method decisions

`docs/adr/0001-differential-expression-and-thresholds.md` records why counts rather than normalised
values go into DESeq2, why genes are reconciled on symbol, why the expression filter is stricter than
the DESeq2 default, and why the log2 fold change threshold is 2 rather than the conventional 1. It
also states the limitation that dataset and biological condition are confounded in this design.

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
committed here, and `data/` is gitignored.

## Layout note

Analysis scripts live in `src/` rather than `scripts/`, which is where the Snakefile expects them.
`scripts/` holds the validation helper only.

## Author

Philipp Graber. [ORCID 0000-0002-3157-1434](https://orcid.org/0000-0002-3157-1434)
