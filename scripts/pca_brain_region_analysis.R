library(tidyverse)
library(readr)
library(readxl)
library(ggpubr)
library(stringr)

# Set file paths for standalone testing
zero_counts_file <- "data/raw/ZERO_RNAseq/ZERO_GeneExpression_rawCounts_CNS_NB_Sarcoma_mRNA_24052024.txt"
normal_counts_file <- "data/raw/PsychENCODE/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.count.txt"
patient_metadata_file <- "data/raw/various/ZERO_psychencode_metadata_new.csv"
patient_metadata_average_file <- "data/raw/various/ZERO_psychencode_metadata_average.xlsx"
housekeeping_genes_file <- "data/raw/Housekeeping_genes/MostStable.csv"

# Output files
pca_all_genes_output <- "results/pca_all_genes.pdf"
pca_housekeeping_output <- "results/pca_housekeeping.pdf"
pca_all_genes_average_output <- "results/pca_all_genes_average.pdf"
pca_housekeeping_average_output <- "results/pca_housekeeping_average.pdf"
boxplots_output_dir <- "results/housekeeping_boxplots"

# Functions from original code
merge_dataframes <- function(ZERO_glioma_counts, normal_counts, join_column_ZERO, join_column_normal) {
  merged_dataframe <- merge(
    ZERO_glioma_counts, 
    normal_counts, 
    by.x = join_column_ZERO, 
    by.y = join_column_normal
  )
  
  # Check merged dataframe
  num_missing_values <- sum(is.na(merged_dataframe))
  cat("Number of missing values in merged dataframe:", num_missing_values, "\n\n")
  
  cat("Dimensions of merged dataframe:", dim(merged_dataframe), "\n\n")
  
  cat("Head of merged dataframe:\n")
  print(head(merged_dataframe))
  
  return(merged_dataframe)
}

prepare_deseq2_analysis <- function(merged_dataframe, patient_metadata) {
  # Convert merged data to matrix and set row names
  DESEQ2_matrix <- as.matrix(merged_dataframe[, -(1:2)])
  rownames(DESEQ2_matrix) <- merged_dataframe$gene_id
  DESEQ2_matrix <- round(DESEQ2_matrix, digits = 0)
  
  # Match patient IDs between DESEQ2 matrix and metadata
  matrix_patient_ids <- colnames(DESEQ2_matrix)
  metadata_patient_ids <- patient_metadata$Patient_ID
  matched_patient_ids <- intersect(matrix_patient_ids, metadata_patient_ids)
  DESEQ2_matrix <- DESEQ2_matrix[, matched_patient_ids]
  
  # Pre-filtering of low count genes
  threshold <- 10  
  keep <- rowSums(DESEQ2_matrix > threshold) >= (0.5 * ncol(DESEQ2_matrix) + 1)
  DESEQ2_matrix <- DESEQ2_matrix[keep, ]
  
  # Reorder metadata to match DESEQ2 matrix column names
  DESEQ2_coldata <- patient_metadata[match(matched_patient_ids, patient_metadata$Patient_ID), ]
  
  # Return list containing DESEQ2 matrix and reordered metadata
  return(list(DESEQ2_matrix = DESEQ2_matrix, DESEQ2_coldata = DESEQ2_coldata))
}

create_pca_plot <- function(data_matrix, col_data, group_var, age_var, pcx = 1, pcy = 2) {
  
  #remove NA
  data_matrix <- na.omit(data_matrix)
  
  # Center and scale the data_matrix for PCA
  scaled_data_matrix <- scale(data_matrix, center = TRUE, scale = TRUE)
  
  # Perform PCA using base R prcomp
  pca_results <- prcomp(t(scaled_data_matrix))
  
  # Extract PCA scores
  pca_scores <- as.data.frame(pca_results$x)
  
  # Add group and age information from col_data
  pca_scores$Group <- col_data[[group_var]]
  pca_scores$Age <- col_data[[age_var]]
  
  # Define colors for the two groups
  my_colors <- c("#1f77b4", "#FF5A3A") 
  
  # Get variance explained by each principal component
  var_explained <- pca_results$sdev^2 / sum(pca_results$sdev^2) * 100
  
  # Create a ggplot
  gg_pca_plot <- ggplot(pca_scores, aes(x = !!sym(paste0("PC", pcx)), y = !!sym(paste0("PC", pcy)), color = Group)) +
    geom_point(shape = 21, size = 3, fill = "white") +  
    scale_color_manual(values = my_colors) +
    labs(x = paste0("PC", pcx, " (", round(var_explained[pcx], 2), "% variance)"), 
         y = paste0("PC", pcy, " (", round(var_explained[pcy], 2), "% variance)"), 
         color = group_var) +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA))
  
  # Return a list containing the ggplot and PCA scores data (matching original structure)
  return(list(gg_pca_plot = gg_pca_plot, pca_scores = pca_scores))
}

# Create output directory for boxplots
if (!dir.exists(boxplots_output_dir)) {
  dir.create(boxplots_output_dir, recursive = TRUE)
}

# Import data (following original approach exactly)
cat("Loading data...\n")
ZERO_counts <- read_tsv(zero_counts_file, show_col_types = FALSE)
normal_counts <- read_tsv(normal_counts_file, show_col_types = FALSE) 
patient_metadata <- read_csv(patient_metadata_file, show_col_types = FALSE)
patient_metadata_average <- read_xlsx(patient_metadata_average_file)
housekeeping_genes <- read_delim(housekeeping_genes_file, delim = ";", escape_double = FALSE, trim_ws = TRUE)
housekeeping_list <- as.list(housekeeping_genes$`Gene name`)

# Process data exactly like original code
normal_counts$Geneid <- sapply(strsplit(normal_counts$Geneid, "\\|"), function(x) x[2])

# Get unique patient IDs from the glioma subset and filter count matrix with unique patient IDs
glioma_patients <- subset(patient_metadata, Type == "Tumor")
unique_patient_ID <- unique(glioma_patients$Patient_ID)

columns_to_keep <- c("gene_id", "transcript_id-s-", unique_patient_ID)
existing_columns <- colnames(ZERO_counts)
valid_columns <- columns_to_keep[columns_to_keep %in% existing_columns]
ZERO_glioma_counts <- ZERO_counts[, valid_columns, drop = FALSE]

# PCA 1: Non-averaged data - Whole gene set
cat("Creating PCA for whole gene set (non-averaged)...\n")
dataframe_whole_gene_set <- merge_dataframes(ZERO_glioma_counts, normal_counts, "gene_id", "Geneid")
analysis_whole_gene_set <- prepare_deseq2_analysis(dataframe_whole_gene_set, patient_metadata)
pca_all_genes <- create_pca_plot(analysis_whole_gene_set$DESEQ2_matrix, analysis_whole_gene_set$DESEQ2_coldata, "Type", "Age")$gg_pca_plot
ggsave(filename = pca_all_genes_output, plot = pca_all_genes, width = 10, height = 8)
cat("Saved PCA all genes (non-averaged) to:", pca_all_genes_output, "\n")

# PCA 2: Non-averaged data - Housekeeping genes only
cat("Creating PCA for housekeeping genes (non-averaged)...\n")
housekeeping_genes_matrix <- as.data.frame(analysis_whole_gene_set$DESEQ2_matrix[rownames(analysis_whole_gene_set$DESEQ2_matrix) %in% housekeeping_list, ])
housekeeping_genes_matrix <- rownames_to_column(housekeeping_genes_matrix, var = "Geneid")

dataframe_housekeeping_genes <- merge_dataframes(ZERO_glioma_counts, housekeeping_genes_matrix, "gene_id", "Geneid")
analysis_housekeeping_genes <- prepare_deseq2_analysis(housekeeping_genes_matrix, patient_metadata)
pca_housekeeping <- create_pca_plot(analysis_housekeeping_genes$DESEQ2_matrix, analysis_housekeeping_genes$DESEQ2_coldata, "Type", "Age")$gg_pca_plot
ggsave(filename = pca_housekeeping_output, plot = pca_housekeeping, width = 10, height = 8)
cat("Saved PCA housekeeping genes (non-averaged) to:", pca_housekeeping_output, "\n")

# Average expression across brain regions (exactly like original)
cat("Averaging expression data across brain regions...\n")
patient_related_columns <- names(normal_counts)[-1]
patient_ids <- sub("\\..*$", "", patient_related_columns)
unique_patient_ids <- unique(patient_ids)

averaged_counts_list <- list(Geneid = normal_counts$Geneid)
for (patient_id in unique_patient_ids) {
    patient_col_indices <- grepl(paste0("^", patient_id, "\\."), names(normal_counts))
    patient_cols <- normal_counts[, patient_col_indices]
    averaged_counts <- rowMeans(patient_cols, na.rm = TRUE)
    averaged_counts_list[[patient_id]] <- averaged_counts
}
averaged_normal_counts <- as.data.frame(averaged_counts_list)

# PCA 3: Averaged data - Whole gene set
cat("Creating PCA for whole gene set (averaged)...\n")
dataframe_whole_gene_set_average <- merge_dataframes(ZERO_glioma_counts, averaged_normal_counts, "gene_id", "Geneid")
analysis_whole_gene_set_average <- prepare_deseq2_analysis(dataframe_whole_gene_set_average, patient_metadata_average)
pca_all_genes_average <- create_pca_plot(analysis_whole_gene_set_average$DESEQ2_matrix, analysis_whole_gene_set_average$DESEQ2_coldata, "Type", "Age")$gg_pca_plot
ggsave(filename = pca_all_genes_average_output, plot = pca_all_genes_average, width = 10, height = 8)
cat("Saved PCA all genes (averaged) to:", pca_all_genes_average_output, "\n")

# PCA 4: Averaged data - Housekeeping genes only
cat("Creating PCA for housekeeping genes (averaged)...\n")
housekeeping_genes_matrix_average <- as.data.frame(analysis_whole_gene_set_average$DESEQ2_matrix[rownames(analysis_whole_gene_set_average$DESEQ2_matrix) %in% housekeeping_list, ])
housekeeping_genes_matrix_average <- rownames_to_column(housekeeping_genes_matrix_average, var = "Geneid")

dataframe_housekeeping_genes_average <- merge_dataframes(ZERO_glioma_counts, housekeeping_genes_matrix_average, "gene_id", "Geneid")
analysis_housekeeping_genes_average <- prepare_deseq2_analysis(housekeeping_genes_matrix_average, patient_metadata_average)
pca_housekeeping_average <- create_pca_plot(analysis_housekeeping_genes_average$DESEQ2_matrix, analysis_housekeeping_genes_average$DESEQ2_coldata, "Type", "Age")$gg_pca_plot
ggsave(filename = pca_housekeeping_average_output, plot = pca_housekeeping_average, width = 10, height = 8)
cat("Saved PCA housekeeping genes (averaged) to:", pca_housekeeping_average_output, "\n")

# Generate boxplots for all housekeeping genes across brain regions
cat("Creating boxplots for housekeeping genes across brain regions...\n")

long_data <- normal_counts %>%
  pivot_longer(-Geneid, names_to = "Sample", values_to = "Expression") %>%
  mutate(Location = str_extract(Sample, "(?<=\\.)\\w+$"))  

housekeeping_data <- long_data %>%
  filter(Geneid %in% housekeeping_list)

plot_count <- 0
for (gene in housekeeping_list) {
  gene_data <- housekeeping_data %>%
    filter(Geneid == gene)
  
  if (nrow(gene_data) > 0) {
    p <- ggplot(gene_data, aes(x = Location, y = Expression, fill = Location)) +
      geom_boxplot() +
      labs(title = paste("Expression of", gene, "across Locations"),
           x = "Location",
           y = "Raw counts") +
      theme_pubr() +
      theme(legend.position = "none",
            strip.background = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8))  
    
    ggsave(filename = file.path(boxplots_output_dir, paste("Boxplot_", gene, ".pdf", sep = "")), 
           plot = p, width = 10, height = 6)
    plot_count <- plot_count + 1
  }
}

cat("Generated", plot_count, "boxplots for housekeeping genes in:", boxplots_output_dir, "\n")
cat("PCA and brain region analysis completed successfully!\n")
