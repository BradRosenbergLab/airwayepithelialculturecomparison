# Returns a seurat object with predictions transfered from a reference dataset
label_transfer_seurat <- function(seurat,
                                  reference,
                                  query_assay,
                                  ref_assay,
                                  norm_method = "pca",
                                  labels = NULL,
                                  name = "ref",
                                  ndims = 30) {
  DefaultAssay(reference) <- ref_assay
  anchors <- FindTransferAnchors(
    reference = reference,
    query = seurat,
    dims = 1:ndims,
    reference.reduction = "pca",
    normalization.method = norm_method,
    query.assay = query_assay,
    reference.assay = ref_assay,
    recompute.residuals = FALSE
  )
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference@meta.data[[labels]],
    dims = 1:ndims
  )
  colnames(predictions) <- paste(name, colnames(predictions), sep = "_")
  seurat <- AddMetaData(seurat, metadata = predictions)
  return(seurat)
}

# Assign major group
# For each cluster, determine the cell type label with the most cells
# Assign this label as the "major group"
assign_major_group <- function(seurat, reference_name, clusters = "seurat_clusters") {
  mat_predictions <- table(
    unlist(seurat[[clusters]]),
    unlist(seurat[[reference_name]])
  ) %>%
    as.matrix()

  assignments <- colnames(mat_predictions)[apply(mat_predictions, 1, which.max)]
  df_assignments <- data.frame(
    cluster = as.character(rownames(mat_predictions)),
    group_major = assignments
  )

  # Add annotations to Seurat obj metadata
  lookup <- c("cluster")
  names(lookup) <- clusters

  df_metadata <- seurat@meta.data %>%
    mutate(clusters = as.character(clusters)) %>%
    tibble::rownames_to_column(var = "cell_barcode") %>%
    left_join(x = ., y = df_assignments, by = lookup)

  rownames(df_metadata) <- df_metadata$cell_barcode

  # Add to Seurat object
  seurat <- AddMetaData(seurat, metadata = df_metadata %>% select(-cell_barcode))
  return(seurat)
}
