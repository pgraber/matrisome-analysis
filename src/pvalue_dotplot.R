library(tidyverse)
library(readr)
library(ggplot2)

# Read the core matrisome DEG results (all significant core matrisome genes)
dataframe_matrisome <- read_csv("output/CoreMatrisome_DESEQ_results_sig.csv", show_col_types = FALSE)

# Define breakpoints for p-value groups
breakpoints <- c(0, 1e-15, 1e-10, 1e-5, Inf)

# Assign each data point to a p-value group
dataframe_matrisome$pvalue_group <- cut(dataframe_matrisome$padj, breaks = breakpoints, labels = c("<1e-15", "1e-15 to 1e-10", "1e-10 to 1e-5", ">1e-5"))

# Define dot sizes for each p-value group (reversed order)
# Define colors for specific matrisome categories
category_colors <- c(
  "ECM Glycoproteins" = "#337EC1",
  "Collagens" = "#E4C538", 
  "Proteoglycans" = "#F3766E"
)

# Create the size plot for all core matrisome DEGs with color coding by category
size_plot <- ggplot(dataframe_matrisome, aes(x = log2FoldChange, y = `Annotated Gene`, 
                                            size = pvalue_group, 
                                            color = `Annotated Matrisome Category`)) +
  geom_point(stroke = 0.5, color = "black") +  
  geom_point(stroke = 0) +  # Colored fill
  scale_size_manual(values = size_mapping) +
  scale_color_manual(values = category_colors) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Log2 Fold Change (Tumor vs Non-tumor)", 
       y = "Genes", 
       size = "Adjusted p-value",
       color = "Matrisome Category",
       title = "Core Matrisome DEGs - P-value Distribution") +
  theme_minimal() +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), labels = c("-4", "-2", "0", "2", "4")) +
  coord_cartesian(xlim = c(-4, 4)) +
  theme(axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

# Display the plot
print(size_plot)

# Save the plot to a file as pdf
ggsave(filename = "output/DEplot_corematrisome_all.pdf", plot = size_plot, width = 10, height = 12)

# Print summary statistics
cat("Summary of Core Matrisome DEGs:\n")
cat("Total genes:", nrow(dataframe_matrisome), "\n")
cat("P-value group distribution:\n")
print(table(dataframe_matrisome$pvalue_group))

# Show top 10 genes by absolute log2FoldChange
cat("\nTop 10 genes by absolute log2FoldChange:\n")
top_genes <- dataframe_matrisome %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(10) %>%
  select(`Annotated Gene`, log2FoldChange, padj, pvalue_group)
print(top_genes)
