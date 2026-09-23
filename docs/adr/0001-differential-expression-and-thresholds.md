# ADR 0001: Differential expression method and significance thresholds

Status: accepted
Date: 2026-09-24 (recorded retrospectively for a workflow first written in 2025)

## Context

The workflow compares bulk RNA-seq from a paediatric tumour cohort against normal brain reference
samples, then asks which of the differentially expressed genes belong to the matrisome. Two method
choices determine every downstream figure, and neither was written down when the code was first
built. This record fixes that.

The two datasets do not share a preparation pipeline. The tumour counts arrive as one column per
sample. The normal brain reference arrives with multiple technical columns per donor and gene
identifiers in a compound `ENSG|SYMBOL` form. Any comparison has to reconcile both before a model is
fitted.

## Decision

**Counts, not normalised values, into DESeq2.** Normal brain technical replicates are averaged per
donor before merging, so each donor contributes one column and donors are not implicitly weighted by
how many times they were sequenced. Averaged values are rounded, because DESeq2's negative binomial
model requires integer counts.

**Gene identifiers are reconciled on symbol.** The normal reference's compound identifier is split
and the symbol taken, because the tumour matrix is symbol-keyed. The merge is an inner join, so only
genes measured in both datasets are tested.

**Expression filter: a gene is kept if it exceeds 10 counts in more than half the samples.** This is
stricter than DESeq2's default independent filtering and is applied before fitting rather than after.
The intent is to remove genes that are absent from one dataset entirely, which would otherwise
produce large fold changes driven by a platform difference rather than by biology.

**Significance: adjusted p at or below 0.05 and an absolute log2 fold change of at least 2.** The
fold-change requirement is deliberate and it is not the conventional 1. A four-fold change was chosen
because the comparison is across two independently generated datasets, so a modest fold change cannot
be separated from a batch effect. The threshold trades sensitivity for the ability to defend each
gene that survives.

## Consequences

The gene set is conservative by construction. Genes with real but moderate differential expression
are not reported, and that is the accepted cost of comparing across datasets that were never designed
to be compared.

Because the thresholds are the whole argument, they are declared once at the top of
`src/differential_expression.R` as `thr` and `lfc_cutoff` rather than being buried at the point of
use, so that changing them is a visible act.

## Alternatives rejected

**Quantile or median-of-ratios normalisation across the two datasets before testing.** Rejected
because it assumes the two datasets are exchangeable, which is exactly the assumption in question.
DESeq2's internal size-factor estimation already handles library-size differences within the model.

**The conventional log2 fold change of 1.** Rejected for the reason above. Across independently
generated datasets it produces a gene list that cannot be defended against a batch-effect objection.

**Batch correction with the dataset as a covariate.** Not possible here. Dataset and biological
condition are perfectly confounded, since every tumour sample comes from one source and every normal
sample from the other. A covariate for dataset would absorb the entire effect of interest. This is
a real limitation of the comparison and it is stated rather than corrected.
