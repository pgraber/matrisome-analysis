renv::load()

library(tidyverse)
library(readr)
library(MatrisomeAnalyzeR)

# Get input and output from snakemake
deseq_file <- snakemake@input[[1]]
output_file_all <- snakemake@output[[1]]
output_file_core <- snakemake@output[[2]]

# Read the DESeq2 results
data <- read_csv(deseq_file, show_col_types = FALSE)

# Annotate genes according to their matrisome status
ECM_annotated <- matriannotate(data = as.data.frame(data), gene.column = "gene_names", species = "human")

# Write all annotated genes
write_csv(ECM_annotated, output_file_all)

# Filter for core matrisome genes only
ECM_core <- ECM_annotated[ECM_annotated$`Annotated Matrisome Division` == "Core matrisome", ]

# Write core matrisome genes
write_csv(ECM_core, output_file_core)