#' Custom DimPlot wrapper to generate figure-ready plots. UMAPs with circle labels.
#'
#' @param object A Seurat object
#' @param group.ident A grouping variable to color points by. If colored by cluster identity, object must be processed by FindClusters.
#' @param colors A vector of colors
#' @param pt.size Size of points on UMAP
#' @param legend.col Number of columns of legend
#' @param label Logical to add labels
#' @param ciclelabels Logical to add circle labels to cluster labels. Works effectively for distinct clusters on UMAP space.
#' @param byrow Legend specification. List variables by row (TRUE) or by column (FALSE).
#' @param inset Logical to specify placement of legend. Inset (TRUE) or out by the side.
#' @param legend.position Placement of the plot. Inherits from inset parameter. Default: top_left
#' @param shape Shape of plotted points
#' @param legend.title Title of legend
#' @param legend.text.align Alignment of legend text
#'
#' @return A publication ready UMAP visualization, where points are colored by group.identity.
#'
#' @examples
#' library(Seurat)
#' data("pbmc_small")
#' CustomUMAP(object = pbmc_small, group.ident = "RNA_snn_res.1", inset = F)
#'
#' @export CustomUMAP

# Alternative, more refined function to add circle labels to UMAP visualizations
CustomUMAP <- function(object,
                       group.ident = "integrated_snn_res.1",
                       colors = legocolors$hex[-1],
                       pt.size = 0.01,
                       legend.col = 4,
                       plot.title = "Unsupervised clustering",
                       ciclelabels = TRUE,
                       byrow = TRUE,
                       label = FALSE,
                       inset = TRUE,
                       shape = 21,
                       legend.title = "Cluster",
                       legend.position = "top left",
                       legend.text.align = 0) {
  library(Seurat)
  library(ggplot2)
  library(legocolors)

  cellnumber <- dim(object)[2]
  group.levels <- levels(object@meta.data[, group.ident])
  p2 <-
    DimPlot(
      object = object,
      pt.size = pt.size,
      label = label,
      group.by = group.ident
    ) +
    scale_color_manual(values = colors, labels = paste0(1:length(group.levels), ": ", group.levels)) +
    scale_fill_manual(values = colors) +
    xlab(expression("UMAP"[1])) +
    ylab(expression("UMAP"[2])) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_text(hjust = 0, vjust = 0),
      axis.title.y = element_text(hjust = 0, vjust = 0),
      legend.title = element_text(face = "bold"),
      legend.text.align = legend.text.align,
      legend.text = element_text(
        color = "black", inherit.blank = T,
        margin = margin(t = 0, r = 0, b = 0, l = -5)
      ),
      legend.key.width = unit(0, "mm"),
      legend.key = element_rect(fill = NA, inherit.blank = T),
      legend.background = element_rect(fill = NA, color = "black", size = 0.25),
      legend.margin = margin(0, 0, 0, 0),
      axis.line = element_blank(),
      panel.border = element_rect(color = "black", fill = NA)
    ) +
    labs(color = legend.title) +
    ggtitle(paste(plot.title, " (", scales::comma_format(digits = 12)(cellnumber), " cells)", sep = "")) +
    guides(colour = guide_legend(
      override.aes = list(
        shape = 21,
        size = 6,
        color = "black",
        fill = colors[1:length(group.levels)]
      ),
      label = T,
      byrow = byrow,
      ncol = legend.col
    ))
  if (ciclelabels) {
    p2 <- custom.LabelClusters(
      plot = p2,
      id = group.ident,
      clusters = levels(object@meta.data[, group.ident]),
      circle.size = 7,
      text.color = "black",
      text.size = 4,
      shape = shape,
      fill = colors[1:length(levels(object@meta.data[, group.ident]))],
      repel = T
    )
  }
  if (inset) {
    p2.gtable <- reposition_legend(
      aplot = p2,
      offset = c(0.001, 0.001),
      position = legend.position
    )

    p2 <- wrap_ggplot_grob(p2.gtable)
  }
  return(p2)
}

# Modify Seurat functions to plot geom_point (shape) and geom_text for label clusters, underlayed with a circles.
# This is key, as it identifies where the labels are plotted on the UMAP
GetXYAesthetics <- function(plot, geom = "GeomPoint", plot.first = TRUE) {
  geoms <- sapply(
    X = plot$layers,
    FUN = function(layer) {
      return(class(x = layer$geom)[1])
    }
  )
  geoms <- which(x = geoms == geom)
  if (length(x = geoms) == 0) {
    stop("Cannot find a geom of class ", geom)
  }
  geoms <- min(geoms)
  if (plot.first) {
    x <- as.character(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)[2]
    y <- as.character(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)[2]
  } else {
    x <- as.character(x = plot$layers[[geoms]]$mapping$x %||% plot$mapping$x)[2]
    y <- as.character(x = plot$layers[[geoms]]$mapping$y %||% plot$mapping$y)[2]
  }
  return(list("x" = x, "y" = y))
}

# Custom function to add circles labels to UMAP visualizations
custom.LabelClusters <- function(
    plot, # Use DimPlot to generate base ggplot to apply function
    id, # The seurat cluster identifier
    clusters = NULL,
    labels = NULL,
    split.by = NULL,
    repel = F,
    colors = colors,
    circle.size = circle.size,
    text.size = text.size,
    ...) {
  xynames <- unlist(x = GetXYAesthetics(plot = plot), use.names = TRUE)
  if (!id %in% colnames(x = plot$data)) {
    stop("Cannot find variable ", id, " in plotting data")
  }
  if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
    warning("Cannot find splitting variable ", id, " in plotting data")
    split.by <- NULL
  }
  data <- plot$data[, c(xynames, id, split.by)]
  possible.clusters <- as.character(x = na.omit(object = unique(x = data[, id])))
  groups <- levels(data[, id]) # clusters %||% as.character(x = na.omit(object = unique(x = data[, id])))
  if (any(!groups %in% possible.clusters)) {
    stop("The following clusters were not found: ", paste(groups[!groups %in% possible.clusters], collapse = ","))
  }
  labels.loc <- lapply(
    X = groups,
    FUN = function(group) {
      data.use <- data[data[, id] == group, , drop = FALSE]
      data.medians <- if (!is.null(x = split.by)) {
        do.call(
          what = "rbind",
          args = lapply(
            X = unique(x = data.use[, split.by]),
            FUN = function(split) {
              medians <- apply(
                X = data.use[data.use[, split.by] == split, xynames, drop = FALSE],
                MARGIN = 2,
                FUN = median,
                na.rm = TRUE
              )
              medians <- as.data.frame(x = t(x = medians))
              medians[, split.by] <- split
              return(medians)
            }
          )
        )
      } else {
        as.data.frame(x = t(x = apply(
          X = data.use[, xynames, drop = FALSE],
          MARGIN = 2,
          FUN = median,
          na.rm = TRUE
        )))
      }
      data.medians[, id] <- factor(group, levels = groups)
      return(data.medians)
    }
  )
  labels.loc <- do.call(what = "rbind", args = labels.loc)
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop("Length of labels (", length(x = labels), ") must be equal to the number of clusters being labeled (", length(x = labels.loc), ").")
  }
  names(x = labels) <- groups
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }
  geom.use <- ifelse(test = repel, yes = geom_point, no = geom_label)
  plot <- plot + geom.use(
    data = labels.loc, size = circle.size, shape = 21, stroke = 0.66, col = "gray17",
    mapping = aes_string(x = xynames["x"], y = xynames["y"]),
    ...
  ) + geom_text(
    size = text.size,
    data = labels.loc, col = "gray17",
    mapping = aes_string(x = xynames["x"], y = xynames["y"], label = labels.loc[[id]] %>% as.numeric()),
    ...
  )
  return(plot)
}
