# matrisome-analysis

A Snakemake workflow for differential expression and extracellular matrix (matrisome) analysis of
bulk RNA-seq, comparing tumour samples against normal brain reference samples.

## Workflow summary

1. **Differential expression.** DESeq2 on tumour versus normal brain counts. Normal-brain technical
   replicates are averaged per donor first, and genes are reconciled on symbol.
2. **Functional enrichment.** g:Profiler results rendered as GOplot chord diagrams for molecular
   function, biological process and cellular component.
3. **Principal component analysis.** PCA across all genes and across housekeeping genes, per sample
   and per averaged brain region, as a check that the tumour and normal cohorts separate on biology
   rather than on library composition.
4. **Matrisome annotation.** Differentially expressed genes annotated with MatrisomeAnalyzeR, split
   into core matrisome and matrisome-associated.
5. **Figures.** Barplots of core matrisome differential expression, and dotplots of effect size and
   significance across matrisome categories.

## Usage

```bash
snakemake -s workflow/Snakefile --cores 4
```

Or in the pinned container, from the repository root:

```bash
docker build -t matrisome-analysis -f env/Dockerfile .
docker run --rm -v "$PWD":/project -w /project matrisome-analysis \
  snakemake -s workflow/Snakefile -c1
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

Written to `results/` and gitignored.

| Output | Contents |
|---|---|
| `DESEQ2_results.csv`, `DESEQ2_results_sig.csv` | All and significant DESeq2 results |
| `Matrisome_DESEQ_results_annotated*.csv` | Results annotated by matrisome division |
| `CoreMatrisome_DESEQ_results_sig.csv` | Significant core matrisome genes |
| `enrichment/GOplot_chord_{MF,BP,CC}.pdf` | Enrichment chord diagrams |
| `CoreMatrisome_DEG_barplot.pdf`, `dotplot_*.pdf`, `pvalue_dotplot_core_matrisome.pdf` | Figures |
| `pca_*.pdf`, `housekeeping_boxplots/` | PCA and housekeeping-gene diagnostics |

## Parameters

Set at the top of `scripts/differential_expression.R`.

| Parameter | Default | Meaning |
|---|---|---|
| `thr` | 10 | A gene is kept if it exceeds this count in more than half the samples |
| `lfc_cutoff` | 2 | Absolute log2 fold change required for significance, alongside adjusted p ≤ 0.05 |

## Environment

R 4.4.1 with Bioconductor 3.19 and Snakemake 9.9.0. All 148 R package versions are pinned in
`renv.lock`, and `env/Dockerfile` builds the exact environment from that lockfile, so a clean
rebuild restores the same versions rather than whatever is current.

## Licence

MIT. See `LICENSE`.

## Author

Philipp Graber. [ORCID 0000-0002-3157-1434](https://orcid.org/0000-0002-3157-1434)
