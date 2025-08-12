renv::load()

library(tidyverse)
library(readr)
library(readxl)
library(DESeq2)
library(here)

data_dir <- here("data")
out_dir  <- here("output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

zero <- read_tsv(file.path(data_dir, "ZERO_RNAseq/ZERO_GeneExpression_rawCounts_CNS_NB_Sarcoma_mRNA_24052024.txt"))
normal <- read_tsv(file.path(data_dir, "PsychENCODE/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.count.txt"))
meta   <- read_xlsx(file.path(data_dir, "various/ZERO_psychencode_metadata_average.xlsx"))

normal$Geneid <- sapply(strsplit(normal$Geneid, "\\|"), function(x) x[2])

# average normal counts by patient (before first ".")
pid_cols <- names(normal)[-1]
pids <- sub("\\..*$", "", pid_cols)
uniq_pids <- unique(pids)
avg_list <- list(Geneid = normal$Geneid)
for (pid in uniq_pids) {
  idx <- grepl(paste0("^", pid, "\\."), names(normal))
  avg_list[[pid]] <- rowMeans(normal[, idx, drop = FALSE], na.rm = TRUE)
}
normal_avg <- as.data.frame(avg_list)

# keep tumor columns from ZERO
tumor_ids <- unique(subset(meta, Type == "Tumor")$Patient_ID)
keep <- intersect(c("gene_id", "transcript_id-s-", tumor_ids), colnames(zero))
zero_tumor <- zero[, keep, drop = FALSE]

# merge tumor + averaged normals
merged <- merge(zero_tumor, normal_avg, by.x = "gene_id", by.y = "Geneid")

# build count matrix + align metadata
counts <- as.matrix(merged[, -(1:2)])
rownames(counts) <- merged$gene_id
counts <- round(counts)

common <- intersect(colnames(counts), meta$Patient_ID)
counts <- counts[, common, drop = FALSE]
coldata <- meta[match(common, meta$Patient_ID), , drop = FALSE]

# prefilter
thr <- 10
keep_genes <- rowSums(counts > thr) >= (0.5 * ncol(counts) + 1)
counts <- counts[keep_genes, , drop = FALSE]

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ Type)
dds <- DESeq(dds)
res <- as.data.frame(results(dds))
res$gene_names <- rownames(res)
res <- res[, c("gene_names", setdiff(names(res), "gene_names"))]
write_csv(res, file.path(out_dir, "DESEQ2_results.csv"))

res_sig <- subset(res, !is.na(padj) & padj <= 0.05 & (log2FoldChange >= 2 | log2FoldChange <= -2))
write_csv(res_sig, file.path(out_dir, "DESEQ2_results_sig.csv"))
