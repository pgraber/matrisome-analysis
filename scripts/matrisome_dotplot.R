renv::load()

library(tidyverse)
library(readr)
library(ggplot2)

# Get input and output from snakemake
matrisome_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

# Read the complete matrisome annotated results (all genes)
DEG_matrisome <- read_csv(matrisome_file, show_col_types = FALSE)

# Count occurrences of each Matrisome Category and arrange by count
category_counts <- DEG_matrisome %>%
  filter(!is.na(`Annotated Matrisome Category`)) %>%
  filter(`Annotated Matrisome Category` != "Non-matrisome") %>%
  group_by(`Annotated Matrisome Category`, `Annotated Matrisome Division`) %>%
  summarise(Count = n(), .groups = "drop") %>%
  arrange(`Annotated Matrisome Division`, (Count)) %>%
  mutate(`Annotated Matrisome Category` = factor(`Annotated Matrisome Category`, levels = unique(.$`Annotated Matrisome Category`)))

# Print the category counts
print("Matrisome category counts:")
print(category_counts)

### Create dotplot for Matrisome categories
matrisome_cleveland_dotplot <- ggplot(category_counts, aes(x = Count, y = `Annotated Matrisome Category`, fill = `Annotated Matrisome Division`)) +
  geom_segment(
    aes(x = 0, xend = Count, y = `Annotated Matrisome Category`, yend = `Annotated Matrisome Category`), 
    color = "black", 
    alpha = 1
  ) +
  geom_dotplot(
    binaxis = "y", 
    stackdir = "center", 
    position = "dodge", 
    dotsize = 1.5, 
    color = "black"
  ) +
  geom_text(
    aes(label = Count), 
    vjust = -1.2,
    hjust = -0.7,
    size = 5, 
    color = "black", 
    show.legend = FALSE
  ) +
  labs(
    title = "Matrisome Genes Only (Complete Dataset)", 
    x = "Counts", 
    y = "Matrisome Category", 
    fill = "Matrisome Division"
  ) + expand_limits(x = c(0, max(category_counts$Count) * 1.35)) +  # Extend x-axis limits
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    legend.position = "right", 
    panel.grid.major.y = element_line(color = "grey", linetype = "dashed"), 
    panel.grid.minor.y = element_blank(), 
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12)
  )

# Display plot 
print(matrisome_cleveland_dotplot)

# Save plot 
ggsave(output_file, plot = matrisome_cleveland_dotplot, 
       width = 12, height = 8, units = "in")

cat("Dotplot saved to:", output_file, "\n")
