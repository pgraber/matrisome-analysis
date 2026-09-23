library(tidyverse)
library(readr)
library(readxl)
library(DESeq2)

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

merged <- merge(zero_tumor, normal_avg, by.x = "gene_id", by.y = "Geneid")

counts <- as.matrix(merged[, -(1:2)])
rownames(counts) <- merged$gene_id
counts <- round(counts)
common <- intersect(colnames(counts), meta$Patient_ID)
counts <- counts[, common, drop = FALSE]
coldata <- meta[match(common, meta$Patient_ID), , drop = FALSE]

keep_genes <- rowSums(counts > thr) >= (0.5 * ncol(counts) + 1)
counts <- counts[keep_genes, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ Type)
dds <- DESeq(dds, quiet = TRUE)
res <- as.data.frame(results(dds))
res$gene_names <- rownames(res)
res <- res[, c("gene_names", setdiff(names(res), "gene_names"))]

readr::write_csv(res, out_all)
res_sig <- subset(res, !is.na(padj) & padj <= 0.05 & (log2FoldChange >= lfc_cutoff | log2FoldChange <= -lfc_cutoff))
readr::write_csv(res_sig, out_sig)

# Create separate up and down regulated files
res_up <- subset(res_sig, log2FoldChange >= lfc_cutoff)
res_down <- subset(res_sig, log2FoldChange <= -lfc_cutoff)

write_csv(res_up, snakemake@output$up)
write_csv(res_down, snakemake@output$down)


