# matrisome-analysis

A Snakemake workflow for differential expression and extracellular matrix (matrisome) analysis of
bulk RNA-seq, comparing tumour samples against normal brain reference samples.

## Workflow summary

1. **Differential expression.** DESeq2 on tumour versus normal brain counts. Normal-brain technical
   replicates are averaged per donor first, and genes are reconciled on symbol.
2. **Functional enrichment.** g:Profiler results rendered as GOplot chord diagrams for molecular
   function, biological process and cellular component.
3. **Matrisome annotation.** Differentially expressed genes annotated with MatrisomeAnalyzeR, split
   into core matrisome and matrisome-associated.
4. **Figures.** Barplots of core matrisome differential expression, and dotplots of effect size and
   significance across matrisome categories.

## Usage

```bash
snakemake --cores 4
```

`rule all` declares the full target set, so any single output can also be requested by name.

## Inputs

Not included in this repository. `data/` is gitignored.

| Input | Format |
|---|---|
| Tumour RNA-seq counts | TSV, genes by samples, from a controlled-access cohort |
| Normal brain RNA-seq counts | TSV, compound `ENSG\|SYMBOL` gene identifiers, multiple columns per donor |
| Sample metadata | XLSX with `Patient_ID` and `Type` columns |
| g:Profiler enrichment exports | CSV, one per ontology |

Paths are declared per rule in the `Snakefile`.

## Outputs

Written to `output/` and gitignored.

| Output | Contents |
|---|---|
| `DESEQ2_results.csv`, `DESEQ2_results_sig.csv` | All and significant DESeq2 results |
| `Matrisome_DESEQ_results_annotated*.csv` | Results annotated by matrisome division |
| `CoreMatrisome_DESEQ_results_sig.csv` | Significant core matrisome genes |
| `enrichment/GOplot_chord_{MF,BP,CC}.pdf` | Enrichment chord diagrams |
| `CoreMatrisome_DEG_barplot.pdf`, `dotplot_*.pdf` | Figures |

## Parameters

Set at the top of `src/differential_expression.R`.

| Parameter | Default | Meaning |
|---|---|---|
| `thr` | 10 | A gene is kept if it exceeds this count in more than half the samples |
| `lfc_cutoff` | 2 | Absolute log2 fold change required for significance, alongside adjusted p ≤ 0.05 |

## Requirements

R with DESeq2, MatrisomeAnalyzeR, GOplot, tidyverse, readxl and readr. Snakemake 7 or later.

## Licence

MIT. See `LICENSE`.

## Author

Philipp Graber. [ORCID 0000-0002-3157-1434](https://orcid.org/0000-0002-3157-1434)
