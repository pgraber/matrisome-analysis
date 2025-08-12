renv::load()

library(tidyverse)
library(circlize)
library(RColorBrewer)

# === Load DESeq2 results ===
DESEQ_results <- read_csv("output/DESEQ2_results.csv")

DESEQ_results_sig <- DESEQ_results %>%
  filter(!is.na(padj) & padj <= 0.05 & abs(log2FoldChange) >= 1)

# Prepare logFC dataframe for colouring
logfc_df <- DESEQ_results %>%
  filter(!is.na(gene_names), !is.na(log2FoldChange)) %>%
  select(gene = gene_names, log2FoldChange)

# === Load g:Profiler results ===
infile  <- "data/gProfiler/gProfiler_hsapiens_12-8-2025_9-23-12 pm__intersections.csv"
outpref <- "output/circos"
dir.create(dirname(outpref), recursive = TRUE, showWarnings = FALSE)

df_raw <- read_csv(infile, show_col_types = FALSE)

if ("adjusted_p_value" %in% names(df_raw)) {
  df_raw$pval <- df_raw$adjusted_p_value
} else {
  df_raw$pval <- NA_real_
}

df <- df_raw %>%
  mutate(term_display = stringr::str_trunc(term_name, 60))

if ("intersections" %in% names(df) && !"intersection" %in% names(df)) {
  df <- rename(df, intersection = intersections)
}

df <- df %>%
  filter(!is.na(intersection), intersection != "")

palette_terms <- function(n) {
  if (n <= 12) RColorBrewer::brewer.pal(max(3, n), "Set3") else
    colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n)
}

# === Build term–gene edges ===
build_edges <- function(dfin, max_genes_per_term = 50, max_total_genes = 100) {
  if (!"pval" %in% names(dfin)) {
    dfin <- mutate(dfin, pval = 1)
  }
  term_gene <- dfin %>%
    arrange(pval, term_display) %>%
    select(term = term_display, intersection) %>%
    separate_rows(intersection, sep = ",") %>%
    mutate(gene = str_trim(intersection)) %>%
    filter(gene != "")
  term_gene_limited <- term_gene %>%
    group_by(term) %>%
    slice_head(n = max_genes_per_term) %>%
    ungroup()
  gene_counts <- table(term_gene_limited$gene)
  if (length(gene_counts) > max_total_genes) {
    top_genes <- names(sort(gene_counts, decreasing = TRUE)[1:max_total_genes])
    term_gene_limited <- filter(term_gene_limited, gene %in% top_genes)
  }
  select(term_gene_limited, term, gene)
}

# === Plotting function with vertical split and log2FC ===
make_chord <- function(
  edges, logfc_df, tag,
  pdf_w = 11, pdf_h = 11, png_w = 3200, png_h = 3200, res = 300,
  transparency = 0.30, rotation_deg = 90  # vertical split
) {
  if (nrow(edges) == 0) return(invisible())

  terms <- unique(edges$term)
  genes <- sort(unique(edges$gene))
  sector_order <- c(terms, genes)

  # Term colours
  term_cols <- setNames(palette_terms(length(terms)), terms)

  # log2FC → colour
  col_fun <- circlize::colorRamp2(
    breaks = c(min(logfc_df$log2FoldChange, na.rm = TRUE),
               0,
               max(logfc_df$log2FoldChange, na.rm = TRUE)),
    colors = c("blue", "white", "red")
  )

  gene_fc <- logfc_df %>%
    filter(gene %in% genes) %>%
    mutate(colour = col_fun(log2FoldChange))

  # Grey for missing genes
  missing_genes <- setdiff(genes, gene_fc$gene)
  if (length(missing_genes) > 0) {
    gene_fc <- bind_rows(
      gene_fc,
      tibble(gene = missing_genes, log2FoldChange = NA, colour = "#BBBBBB")
    )
  }

  gene_cols <- setNames(gene_fc$colour, gene_fc$gene)
  grid_cols <- c(term_cols, gene_cols)

  # Gap between blocks
  gap_terms <- rep(2, length(terms))
  if (length(gap_terms)) gap_terms[length(gap_terms)] <- 8
  gaps <- c(gap_terms, rep(1, length(genes)))

  plot_fun <- function() {
    circos.clear()
    circos.par(start.degree = rotation_deg, clock.wise = TRUE, gap.after = gaps)
    chordDiagram(
      x = edges,
      order = sector_order,
      grid.col = grid_cols,
      transparency = transparency,
      annotationTrack = NULL,
      preAllocateTracks = list(track.height = 0.06)
    )
    circos.track(track.index = 1, panel.fun = function(x, y) {
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      nm <- get.cell.meta.data("sector.index")
      if (nm %in% genes) {
        circos.text(mean(xlim), ylim[1], nm, facing = "clockwise",
                    niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
      } else {
        circos.text(mean(xlim), ylim[1], nm, facing = "inside",
                    niceFacing = TRUE, adj = c(0.5, 0), cex = 0.9)
      }
    }, bg.border = NA)
  }

  # PDF
  pdf(file.path(dirname(outpref), paste0(basename(outpref), "_", tag, ".pdf")), pdf_w, pdf_h)
  plot_fun()
  dev.off()

  # PNG
  png(file.path(dirname(outpref), paste0(basename(outpref), "_", tag, ".png")), png_w, png_h, res = res)
  plot_fun()
  dev.off()
}

# === Ontology mapping and plotting ===
ont_map <- c(MF = "GO:MF", BP = "GO:BP", CC = "GO:CC")

for (nm in names(ont_map)) {
  cat_code <- ont_map[[nm]]
  cat_df <- filter(df, source == cat_code)
  cat_edges <- build_edges(cat_df, max_genes_per_term = 15, max_total_genes = 60)
  n_terms <- length(unique(cat_edges$term))
  n_genes <- length(unique(cat_edges$gene))
  message("Ontology ", nm, ": terms=", n_terms, " genes=", n_genes, " edges=", nrow(cat_edges))
  make_chord(cat_edges, logfc_df, nm)
}
