renv::load()

library(tidyverse)
library(readr)
library(ggplot2)
library(ggpubr)

# Get input and output from snakemake
deseq_sig_file <- snakemake@input$deseq_sig
matrisome_file <- snakemake@input$matrisome
output_file <- snakemake@output[[1]]

# Read the DESeq2 results and matrisome annotated results
DESEQ_results_sig <- read_csv(deseq_sig_file, show_col_types = FALSE)
DEG_core_matrisome <- read_csv(matrisome_file, show_col_types = FALSE)

# Create stacked bar chart to indicate DEG results 

# Define counts
total_DEG <- length(unique(DESEQ_results_sig$gene_names))
matrisome_DEG <- length(unique(DEG_core_matrisome$`Annotated Gene`))

# Define data frame
data <- data.frame(
  Category = c("Total DEG", "Core matrisome DEG"),
  Count = c(total_DEG, matrisome_DEG)
)

# Create stacked bar plot
(barplot_DEG <- ggplot(data, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.2) +
  coord_flip() +
  labs(title = paste("Counts of Unique Genes:", paste(data$Count, collapse = " / ")),
       x = NULL,
       y = NULL) + # Remove y-axis label
  theme_pubr() +
  scale_fill_manual(values = c("Total DEG" = "blue", "Core matrisome DEG" = "red")) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "top"))

# Save the plot
ggsave(plot = barplot_DEG, filename = output_file, 
       width = 8, height = 6, units = "in")

