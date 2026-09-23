# Load environment and libraries
renv::load()

library(tidyverse)
library(circlize)
library(RColorBrewer)
library(GOplot)
library(viridis)
library(here)

# Input files from Snakemake
deseq_file <- snakemake@input[["deseq_results"]]
gprofiler_file <- snakemake@input[["gprofiler_data"]]

# Output files from Snakemake
chord_output <- snakemake@output[["chord_plot"]]
circ_output <- snakemake@output[["circ_data"]]

# Get ontology from wildcards for title
ontology <- snakemake@wildcards[["ontology"]]

# Create output directory
output_dir <- dirname(chord_output)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load DESeq2 results
DESEQ_results <- read_csv(deseq_file)

DESEQ_results_sig <- DESEQ_results %>%
  filter(!is.na(padj) & padj <= 0.05 & abs(log2FoldChange) >= 1)

# Read all gProfiler files to calculate global limits from combined datasets
gprofiler_files <- c(
  "data/raw/gProfiler/gProfiler_hsapiens_26-8-2025_8-56-28 am__intersections_MF.csv",
  "data/raw/gProfiler/gProfiler_hsapiens_26-8-2025_8-56-28 am__intersections_BP.csv",
  "data/raw/gProfiler/gProfiler_hsapiens_26-8-2025_8-56-28 am__intersections.csv"
)

# Function to process each gProfiler file and get genes with logFC
process_gprofiler_file <- function(file_path) {
  raw_data <- read_csv(file_path)
  
  # Select and rename columns
  go_data <- raw_data %>%
    select(
      Category = source,
      ID = term_id,
      term = term_name,
      Genes = intersections,
      adj_pval = adjusted_p_value
    )
  
  # Expand gene lists
  go_expanded <- go_data %>%
    separate_rows(Genes, sep = ",") %>%
    filter(!is.na(Genes) & Genes != "")
  
  # Join with significant genes
  circ <- go_expanded %>%
    inner_join(DESEQ_results_sig, by = c("Genes" = "gene_names")) %>%
    select(
      Category,
      ID,
      term,
      genes = Genes,
      adj_pval,
      logFC = log2FoldChange
    )
  
  return(circ)
}

# Process all files and combine
all_go_data <- map_dfr(gprofiler_files, process_gprofiler_file)

# Calculate global limits from all GO genes across all ontologies
global_min <- min(all_go_data$logFC, na.rm = TRUE)
global_max <- max(all_go_data$logFC, na.rm = TRUE)
global_limits <- c(global_min, global_max)

# Load g:Profiler results for current ontology
raw_data <- read_csv(gprofiler_file)

### Prepare data for GOPlot input ###

# Select and rename columns
go_data <- raw_data %>%
  select(
    Category = source,
    ID = term_id,
    term = term_name,
    Genes = intersections,
    adj_pval = adjusted_p_value
  )

# Expand gene lists
go_expanded <- go_data %>%
  separate_longer_delim(Genes, delim = ",") %>%
  mutate(Genes = str_trim(Genes)) %>%
  filter(!is.na(Genes), Genes != "")

# Join with all significant genes
circ <- go_expanded %>%
  left_join(
    DESEQ_results_sig %>% 
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

# Save circ object
write_csv(circ, circ_output)

###

# Filter to top genes per term
circ_filtered <- circ %>%
  group_by(term) %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 12) %>%
  ungroup()

# Create chord matrix
process_list <- unique(circ_filtered$term)
chord <- chord_dat(data = circ_filtered, process = process_list)

# Add logFC values to matrix
gene_logfc <- circ_filtered %>%
  select(genes, logFC) %>%
  distinct()

chord_final <- chord %>%
  as.data.frame() %>%
  rownames_to_column("genes") %>%
  left_join(gene_logfc, by = "genes") %>%
  column_to_rownames("genes")

# Create and save chord plot
chord_plot <- GOChord(
  data = chord_final,
  title = paste0('GO ', ontology, ' Enrichment (Top 15 genes per term)'),
  space = 0.01,
  gene.order = 'logFC',
  gene.space = 0.2,
  gene.size = 6,
  nlfc = 1,
  lfc.col = c("blue", "white", "red"),
  lfc.min = 0,
  lfc.max = 100,
  ribbon.col = viridis::viridis(length(process_list), option = "D"),
  border.size = 0,    
  process.label = 8
)

# Add your custom red color scale with global limits
chord_plot_custom <- chord_plot +
  scale_fill_gradient(
    low = "#FFE6E6",    
    high = "red",
    name = "log2FC",
    limits = global_limits
  ) 

print(chord_plot_custom)

# Save output
ggsave(
  filename = chord_output,
  plot = chord_plot_custom,
  width = 18,
  height = 18
)