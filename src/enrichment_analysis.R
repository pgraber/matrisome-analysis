renv::load()

library(tidyverse)
library(circlize)
library(RColorBrewer)
library(GOplot)

# === Load DESeq2 results ===
DESEQ_results <- read_csv("output/DESEQ2_results.csv")

DESEQ_results_sig <- DESEQ_results %>%
  filter(!is.na(padj) & padj <= 0.05 & abs(log2FoldChange) >= 1)

# === Load g:Profiler results ===
infile  <- "data/gProfiler/gProfiler_hsapiens_12-8-2025_9-23-12 pm__intersections_MF.csv"
outpref <- "output/circos"
dir.create(dirname(outpref), recursive = TRUE, showWarnings = FALSE)

raw_data <- read_csv(infile)

# === Prepare data for GOplot ===
circ <- raw_data %>%
  select(
    Category = source,
    ID = term_id, 
    term = term_name,
    Genes = intersections, 
    adj_pval = adjusted_p_value
  ) %>%
  separate_longer_delim(Genes, delim = ",") %>%
  mutate(Genes = str_trim(Genes)) %>%
  filter(!is.na(Genes), Genes != "") %>%
  # Join with upregulated genes only
  left_join(
    DESEQ_results_sig %>% 
      filter(log2FoldChange > 0) %>%
      select(gene_names, log2FoldChange),
    by = c("Genes" = "gene_names")
  ) %>%
  filter(!is.na(log2FoldChange)) %>%
  select(
    category = Category,
    ID = ID,
    term = term,
    genes = Genes,
    adj_pval = adj_pval,
    logFC = log2FoldChange
  )

print("Original circ object:")
print(paste("Total gene-term pairs:", nrow(circ)))

# === Filter to reduce number of genes ===
circ_filtered <- circ %>%
  group_by(term) %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 15) %>%  # Top 15 genes per term
  ungroup()

print("Filtered circ object:")
print(paste("Reduced to gene-term pairs:", nrow(circ_filtered)))
print(paste("Unique genes:", length(unique(circ_filtered$genes))))

# === Create chord matrix with filtered data ===
process_list <- unique(circ_filtered$term)
chord <- chord_dat(data = circ_filtered, process = process_list)

# Add logFC values to chord matrix
gene_logfc <- circ_filtered %>%
  select(genes, logFC) %>%
  distinct()

chord_final <- chord %>%
  as.data.frame() %>%
  rownames_to_column("genes") %>%
  left_join(gene_logfc, by = "genes") %>%
  column_to_rownames("genes")

print("Final chord matrix:")
print(paste("Matrix dimensions:", nrow(chord_final), "x", ncol(chord_final)))

# === Create and save chord plot ===
chord_plot <- GOChord(
  data = chord_final, 
  title = 'GO Molecular Function Enrichment (Top 15 genes per term)',
  space = 0.02, 
  gene.order = 'logFC',
  gene.space = 0.25, 
  gene.size = 4,
  nlfc = 1,
  ribbon.col = brewer.pal(length(process_list), "Set2"),
  border.size = 0.5,
  process.label = 8
)

print(chord_plot)

ggsave(
  filename = paste0(outpref, "_GOplot_chord_filtered.pdf"),
  plot = chord_plot,
  width = 12,
  height = 12
)

print("Filtered chord plot created and saved!")