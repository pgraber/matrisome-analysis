library(tidyverse)
library(readr)
library(MatrisomeAnalyzeR)

# Get input and output from snakemake
deseq_all_file <- snakemake@input$all_results
deseq_sig_file <- snakemake@input$sig_results
output_all_annotated <- snakemake@output$all_annotated
output_sig_annotated <- snakemake@output$sig_annotated
output_core_sig <- snakemake@output$core_sig

# Read the full DESeq2 results
data_all <- read_csv(deseq_all_file, show_col_types = FALSE)

# Annotate all genes according to their matrisome status
ECM_all_annotated <- matriannotate(data = as.data.frame(data_all), gene.column = "gene_names", species = "human")

# Write all annotated genes (complete dataset)
write_csv(ECM_all_annotated, output_all_annotated)

# Read the significant DESeq2 results
data_sig <- read_csv(deseq_sig_file, show_col_types = FALSE)

# Annotate significant genes according to their matrisome status
ECM_sig_annotated <- matriannotate(data = as.data.frame(data_sig), gene.column = "gene_names", species = "human")

# Write all significant annotated genes
write_csv(ECM_sig_annotated, output_sig_annotated)

# Filter for core matrisome genes only (from significant results)
ECM_core <- ECM_sig_annotated[ECM_sig_annotated$`Annotated Matrisome Division` == "Core matrisome", ]

# Write core matrisome genes (from significant results)
write_csv(ECM_core, output_core_sig)