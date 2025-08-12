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
# Expand genes and join with expression data
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

print("Circ object ready:")
head(circ)

# === Create chord matrix ===
process_list <- unique(circ$term)
chord <- chord_dat(data = circ, process = process_list)

# Add logFC values to chord matrix
gene_logfc <- circ %>%
  select(genes, logFC) %>%
  distinct()

chord_final <- chord %>%
  as.data.frame() %>%
  rownames_to_column("genes") %>%
  left_join(gene_logfc, by = "genes") %>%
  column_to_rownames("genes")

print("Final chord matrix:")
head(chord_final)

# === Create and save chord plot ===
chord_plot <- GOChord(
  data = chord_final, 
  title = 'GO Molecular Function Enrichment',
  space = 0.02, 
  gene.order = 'logFC',
  gene.space = 0.25, 
  gene.size = 4,
  nlfc = 1,
  ribbon.col = brewer.pal(length(process_list), "Set2"),
  border.size = 0.5,
  process.label = 8, 
  limit = c(20,0)
)

print(chord_plot)

ggsave(
  filename = paste0(outpref, "_GOplot_chord.png"),
  plot = chord_plot,
  width = 14,
  height = 12,
  dpi = 300,
  bg = "white"
)

print("Chord plot created and saved!")