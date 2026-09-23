library(tidyverse)
library(readr)
library(readxl)
library(DESeq2)

source("scripts/validate.R")
vcheck_init("output/checks", step = "differential_expression")

zero_file   <- snakemake@input$zero_counts
normal_file <- snakemake@input$normal_counts
meta_file   <- snakemake@input$metadata
out_all     <- snakemake@output$all
out_sig     <- snakemake@output$sig
thr         <- 10
lfc_cutoff  <- 2

zero   <- read_tsv(zero_file, show_col_types = FALSE)
normal <- read_tsv(normal_file, show_col_types = FALSE)
meta   <- read_xlsx(meta_file)

normal$Geneid <- sapply(strsplit(normal$Geneid, "\\|"), function(x) if (length(x) >= 2) x[2] else x[1])

pid_cols <- names(normal)[-1]
pids <- sub("\\..*$", "", pid_cols)
uniq_pids <- unique(pids)
avg_list <- list(Geneid = normal$Geneid)
for (pid in uniq_pids) {
  idx <- grepl(paste0("^", pid, "\\."), names(normal))
  avg_list[[pid]] <- rowMeans(normal[, idx, drop = FALSE], na.rm = TRUE)
}
normal_avg <- as.data.frame(avg_list)

tumor_ids <- unique(subset(meta, Type == "Tumor")$Patient_ID)
keep <- intersect(c("gene_id", "transcript_id-s-", tumor_ids), colnames(zero))
zero_tumor <- zero[, keep, drop = FALSE]

vcheck("tumour samples selected", length(intersect(tumor_ids, colnames(zero))),
       length(intersect(unique(subset(meta, Type == "Tumor")$Patient_ID), colnames(zero))),
       source = "metadata sheet, Type == Tumor")

genes_expected <- length(intersect(zero_tumor$gene_id, normal_avg$Geneid))
merged <- merge(zero_tumor, normal_avg, by.x = "gene_id", by.y = "Geneid")
vcheck("genes retained after join", nrow(merged), genes_expected,
       source = "intersection of tumour and normal gene IDs, computed before the join")
vcheck("no duplicate gene IDs after join", !any(duplicated(merged$gene_id)), TRUE,
       source = "one row per gene by design")

counts <- as.matrix(merged[, -(1:2)])
rownames(counts) <- merged$gene_id
counts <- round(counts)
common <- intersect(colnames(counts), meta$Patient_ID)
counts <- counts[, common, drop = FALSE]
coldata <- meta[match(common, meta$Patient_ID), , drop = FALSE]

keep_genes <- rowSums(counts > thr) >= (0.5 * ncol(counts) + 1)
counts <- counts[keep_genes, , drop = FALSE]

vcheck("samples matched to metadata", nrow(coldata), ncol(counts),
       source = "count matrix columns")
vcheck("no negative counts", min(counts) >= 0, TRUE,
       source = "counts are non-negative by definition")

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ Type)
dds <- DESeq(dds, quiet = TRUE)
res <- as.data.frame(results(dds))
res$gene_names <- rownames(res)
res <- res[, c("gene_names", setdiff(names(res), "gene_names"))]

vcheck("results rows equal genes tested", nrow(res), nrow(counts),
       source = "genes surviving the expression filter")
vcheck("p-values within [0, 1]",
       all(is.na(res$pvalue) | (res$pvalue >= 0 & res$pvalue <= 1)), TRUE,
       source = "definition of a p-value")

readr::write_csv(res, out_all)
res_sig <- subset(res, !is.na(padj) & padj <= 0.05 & (log2FoldChange >= lfc_cutoff | log2FoldChange <= -lfc_cutoff))
vcheck("significant genes are a subset of all genes", nrow(res_sig) <= nrow(res), TRUE,
       source = "arithmetic")
readr::write_csv(res_sig, out_sig)
vcheck_summary()
