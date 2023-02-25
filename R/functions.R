require(scDblFinder)
require(scuttle)
require(Seurat)
require(ggplot2)
require(ggrepel)
require(pals)
require(scales)

# Returns a set of n colors. Identical to ggplot2 defaults.
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Wrapper function for running scDblFinder workflow. Returns a seurat object.
detect_doublets <- function(seurat, clusters = NULL, knownDoublets = NULL) {
  # Convert to sce object for scDblFinder workflow
  DefaultAssay(seurat) <- "RNA"
  ##  Add metadata column for HTO classification (logical)
  seurat$Doublet <- seurat$HTO_classification.global == "Doublet"
  sce <- as.SingleCellExperiment(seurat, assay = "RNA")

  # Run scDblFinder
  ## Primary use defaults
  sce <- scDblFinder(
    sce = sce,
    clusters = NULL, # Generate artificial doublets
    knownDoublets = knownDoublets # Specify HTO classifications
  )
  # Add scDblFinder results as metadata back to Seurat object
  df_metadata <- as.data.frame(
    colData(sce)[, c(
      "ident", "scDblFinder.class", "scDblFinder.score",
      "scDblFinder.weighted", "scDblFinder.cxds_score"
    )]
  )
  seurat <- AddMetaData(seurat, metadata = df_metadata)
  return(seurat)
}

# Wrapper for running tSNE embedding on HTO assay distances
run_tSNE_HTO <- function(seurat, remove.negatives = TRUE, perplexity = 100) {
  Idents(seurat) <- "hash.ID"
  hto_tsne <- seurat

  # Remove droplets designated "Negative" by demultiplexing?
  if (remove.negatives == TRUE) {
    hto_tsne <- subset(hto_tsne, idents = "Negative", invert = TRUE)
  }

  # Calculate a distance matrix using HTO
  hto.dist.mtx <- as.matrix(
    dist(t(GetAssayData(object = hto_tsne, assay = "HTO")))
  )

  # Calculate tSNE embeddings with a distance matrix
  hto_tsne <- RunTSNE(hto_tsne,
    distance.matrix = hto.dist.mtx,
    perplexity = perplexity
  )
  return(hto_tsne)
}

# Wrapper for drawing a PCA scatterplot with ggplot2
drawPCAPlot <- function(
    pca.dat,
    rel_variances = rel.var.mat,
    x = "PC1",
    y = "PC2",
    color = "ident",
    shape = "SampleGroup") {
  ggplot(
    pca.dat,
    aes_string(x, y, color = color, shape = shape),
  ) +
    geom_point(size = 3) +
    theme_classic() +
    scale_color_hue(drop = FALSE) +
    labs(
      x = paste0(x, " (", round(rel_variances[x] * 100, digits = 1), "%)"),
      y = paste0(y, " (", round(rel_variances[y] * 100, digits = 1), "%)")
    )
}

# A 'primitive' volcano plot. Expects columns `logFC` and `PValue` in `data`.
# Labels will be annotated for all genes in `to.label` vector. Use cmap to pass
# color mapping, will default to standard ggplot colors.
myVolcano <- function(data, to.label = NULL, cmap = NULL, alpha_min = 0.01, label_size = 4, ...) {
  if (is.null(cmap)) {
    cmap <- gg_color_hue(data %>%
      {
        n_distinct(.$color)
      })
  }
  p <- ggplot(
    data,
    aes(
      logFC,
      -log10(PValue),
      ...
    )
  ) +
    geom_point(size = 1) +
    geom_text_repel(data = data %>% filter({{ to.label }}), aes(label = gene), color = "black", box.padding = 0.5, size = label_size) +
    scale_color_manual(limits = names(cmap), values = cmap) +
    scale_alpha(range = c(alpha_min, 1)) +
    theme_bw()
  return(p)
}

# A 'decorated' volcano plot highlighting significant genes in both directions.
# Labels the top n genes in each direction with `label.n`. Option to pass
# additional set labels in `manual_ids`. Expects columns "PValue", "logFC" and
# "is.deg" to exist in `data`.
UpDownVolcano <- function(data, label.n = 10, alpha = 0.6, manual_ids = NULL) {
  # select top n points
  data <- data %>%
    arrange(PValue) %>%
    group_by(logFC < 0) %>%
    mutate(
      n = 1:n(),
      to.label = case_when(
        is.deg == 0 ~ FALSE,
        n <= label.n ~ TRUE,
        TRUE ~ FALSE
      ),
      to.label = if_else(gene %in% manual_ids, TRUE, to.label)
    )
  # plot call
  myVolcano(
    data,
    color = factor(.data[["is.deg"]], ),
    to.label = to.label,
    cmap = c("-1" = "blue", "1" = "red", "0" = "grey60"),
    alpha = I(alpha)
  ) +
    guides(color = "none", alpha = "none")
}

# Highlight a set of genes in a volcano. Expects "PValue" and "logFC" to exist in
# `de.results`. The "alpha_min" argument helps accentuate the selected set.
highlightVolcanoByGene <- function(de.results, ids, alpha_min = 0.01) {
  de.results.filt <- de.results %>%
    mutate(
      to.label = gene %in% ids,
      color = ifelse(to.label, "yes", "no"),
    ) %>%
    arrange(color)

  p <- myVolcano(
    de.results.filt,
    to.label = to.label,
    cmap = c("yes" = "red", "no" = "grey60"),
    color = .data[["color"]],
    alpha = .data[["to.label"]] %>% as.integer(),
    alpha_min = alpha_min
  ) +
    guides(color = "none", alpha = "none")
  p
}

# Scalable ggplot theme for dot plot
dot_theme_fxn <- function(base_size = 14) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      # Figure title
      plot.title = element_text(
        size = rel(1), face = "bold",
        margin = margin(0, 0, 5, 0), hjust = 0
      ),
      # Plotting area
      panel.border = element_rect(size = 0.25, fill = "transparent"),
      panel.grid = element_blank(),
      # Axis labels
      axis.title = element_blank(),
      axis.text = element_text(size = rel(0.60), face = "plain"),
      axis.text.x = element_text(angle = 60, hjust = 0, vjust = 0),
      axis.line = element_blank(), # axis lines from panel border
      # axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(size = 0.25),
      # Legend
      legend.title = element_text(size = rel(0.66), face = "plain"),
      legend.text = element_text(size = rel(0.66), face = "plain"),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(5, "mm"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin()
    )
}

# A simple edgeR fit wrapper
run_standard_edgeR_fit <- function(counts, group, design) {
  y <- DGEList(counts = counts)
  # Filter by expression
  keep <- filterByExpr(y, group = group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  # Normalize
  y <- calcNormFactors(y)
  # Dispersion
  y <- estimateDisp(y, design)
  # Fit
  fit <- glmQLFit(y, design)
  return(fit)
}
