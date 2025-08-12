renv::load()

library(tidyverse)   # includes dplyr, readr, stringr, etc.
library(circlize)
library(RColorBrewer)

# Removed clusterProfiler / org.Hs.eg.db / ReactomePA / here to avoid masking dplyr::select
# (AnnotationDbi::select was capturing tibble input and causing the error).

DESEQ_results <- read_csv("output/DESEQ2_results.csv")

DESEQ_results_sig <- DESEQ_results %>%
    filter(!is.na(padj) & padj <= 0.05 & abs(log2FoldChange) >= 1)
DESEQ_up <- DESEQ_results_sig %>% filter(log2FoldChange > 0)
DESEQ_down <- DESEQ_results_sig %>% filter(log2FoldChange < 0)

write_csv(DESEQ_up, "output/DESEQ2_results_up.csv")
write_csv(DESEQ_down, "output/DESEQ2_results_down.csv")


infile  <- "data/gProfiler/gProfiler_hsapiens_12-8-2025_9-23-12 pm__intersections.csv"
outpref <- "output/circos"
dir.create(dirname(outpref), recursive = TRUE, showWarnings = FALSE)

df_raw <- readr::read_csv(infile, show_col_types = FALSE)

# Use adjusted_p_value column directly (user confirmed this is present)
if ("adjusted_p_value" %in% names(df_raw)) {
  df_raw$pval <- df_raw$adjusted_p_value
} else {
  df_raw$pval <- NA_real_
}

df <- df_raw %>%
  mutate(term_display = stringr::str_trunc(term_name, 60))

# Harmonize intersections column name
if ("intersections" %in% names(df) && !"intersection" %in% names(df)) {
  df <- dplyr::rename(df, intersection = intersections)
}

df <- df %>%
  filter(!is.na(intersection), intersection != "")

palette_terms <- function(n) {
  if (n <= 12) RColorBrewer::brewer.pal(max(3, n), "Set3") else
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n)
}

# Build term-term adjacency matrix (shared gene counts) for summary chord
build_term_adjacency <- function(dfin) {
  if (nrow(dfin) == 0) return(NULL)
  tg <- dfin %>%
    dplyr::select(term_display, intersection) %>%
    tidyr::separate_rows(intersection, sep = ",") %>%
    dplyr::mutate(gene = stringr::str_trim(intersection)) %>%
    dplyr::filter(gene != "") %>%
    dplyr::distinct(term_display, gene)
  terms <- unique(tg$term_display)
  if (length(terms) <= 1) return(NULL)
  term_genes <- split(tg$gene, factor(tg$term_display, levels = terms))  # ensure order
  m <- matrix(0, nrow = length(terms), ncol = length(terms), dimnames = list(terms, terms))
  n <- length(terms)
  if (n < 2) return(NULL)
  for (i in seq_len(n - 1)) {
    gi <- term_genes[[i]]
    for (j in (i + 1):n) {
      gj <- term_genes[[j]]
      shared <- length(intersect(gi, gj))
      m[i, j] <- shared
      m[j, i] <- shared
    }
  }
  # Remove rows/cols that are all zero (no sharing)
  keep <- which(rowSums(m) > 0)
  if (length(keep) < 2) return(NULL)
  m[keep, keep, drop = FALSE]
}

build_edges <- function(dfin, max_genes_per_term = 50, max_total_genes = 100) {
  # ensure pval exists
  if (!"pval" %in% names(dfin)) {
    dfin <- dfin %>% mutate(pval = 1)
  }
  
  # Create term-gene pairs
  term_gene <- dfin %>%
    dplyr::arrange(pval, term_display) %>%
    dplyr::select(term = term_display, intersection) %>%
    tidyr::separate_rows(intersection, sep = ",") %>%
    dplyr::mutate(gene = stringr::str_trim(intersection)) %>%
    dplyr::filter(gene != "")
  
  # Limit genes per term to prevent overcrowding
  term_gene_limited <- term_gene %>%
    dplyr::group_by(term) %>%
    dplyr::slice_head(n = max_genes_per_term) %>%
    dplyr::ungroup()
  
  # If still too many unique genes, keep only most frequent ones
  gene_counts <- table(term_gene_limited$gene)
  if (length(gene_counts) > max_total_genes) {
    top_genes <- names(sort(gene_counts, decreasing = TRUE)[1:max_total_genes])
    term_gene_limited <- term_gene_limited %>%
      dplyr::filter(gene %in% top_genes)
  }
  
  # Return all term-gene pairs (bipartite structure)  
  term_gene_limited %>%
    dplyr::select(term, gene)
}

make_chord <- function(edges, tag, pdf_w = 11, pdf_h = 11, png_w = 3200, png_h = 3200, res = 300, transparency = 0.30) {
  if (nrow(edges) == 0) return(invisible())
  terms <- unique(edges$term)
  genes <- sort(unique(edges$gene))
  grid_cols <- c(
    setNames(palette_terms(length(terms)), terms),
    setNames(rep("#BBBBBB", length(genes)), genes)
  )
  
  # PDF version
  pdf(file.path(dirname(outpref), paste0(basename(outpref), "_", tag, ".pdf")), pdf_w, pdf_h)
  circos.clear()
  chordDiagram(edges,
               grid.col = grid_cols,
               transparency = transparency,
               annotationTrack = NULL,  # Remove default labels
               preAllocateTracks = list(track.height = 0.06))
  # Add custom rotated labels for genes (pointing outward)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    # Rotate gene labels 90 degrees outward
    if (sector.name %in% genes) {
      circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", 
                  niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
    } else {
      # Keep term labels horizontal
      circos.text(mean(xlim), ylim[1], sector.name, facing = "inside", 
                  niceFacing = TRUE, adj = c(0.5, 0), cex = 0.9)
    }
  }, bg.border = NA)
  dev.off()
  
  # PNG version
  png(file.path(dirname(outpref), paste0(basename(outpref), "_", tag, ".png")), png_w, png_h, res = res)
  circos.clear()
  chordDiagram(edges,
               grid.col = grid_cols,
               transparency = transparency,
               annotationTrack = NULL,  # Remove default labels
               preAllocateTracks = list(track.height = 0.06))
  # Add custom rotated labels for genes (pointing outward)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    # Rotate gene labels 90 degrees outward
    if (sector.name %in% genes) {
      circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", 
                  niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
    } else {
      # Keep term labels horizontal
      circos.text(mean(xlim), ylim[1], sector.name, facing = "inside", 
                  niceFacing = TRUE, adj = c(0.5, 0), cex = 0.9)
    }
  }, bg.border = NA)
  dev.off()
}

make_term_chord <- function(adj, tag, pdf_w = 8, pdf_h = 8, png_w = 2000, png_h = 2000, res = 300, transparency = 0.25) {
  if (is.null(adj)) return(invisible())
  terms <- rownames(adj)
  grid_cols <- setNames(palette_terms(length(terms)), terms)
  pdf(file.path(dirname(outpref), paste0(basename(outpref), "_", tag, "_summary.pdf")), pdf_w, pdf_h)
  circos.clear()
  chordDiagram(adj,
               grid.col = grid_cols,
               transparency = transparency,
               annotationTrack = "name",
               directional = FALSE)
  dev.off()
  png(file.path(dirname(outpref), paste0(basename(outpref), "_", tag, "_summary.png")), png_w, png_h, res = res)
  circos.clear()
  chordDiagram(adj,
               grid.col = grid_cols,
               transparency = transparency,
               annotationTrack = "name",
               directional = FALSE)
  dev.off()
}

ont_map <- c(MF = "GO:MF", BP = "GO:BP", CC = "GO:CC")

for (nm in names(ont_map)) {
  cat_code <- ont_map[[nm]]
  cat_df <- dplyr::filter(df, source == cat_code)
  cat_edges <- build_edges(cat_df, max_genes_per_term = 25, max_total_genes = 70)
  n_terms <- length(unique(cat_edges$term))
  n_genes <- length(unique(cat_edges$gene))
  message("Ontology ", nm, ": terms=", n_terms, " genes=", n_genes, " edges=", nrow(cat_edges))
  # Always create term-gene bipartite chord diagram
  make_chord(cat_edges, nm)
}

