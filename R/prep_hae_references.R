# prep_hae_references.R
# Read in and prepare reference scRNA-Seq HAE datasets
# If run for the first time on minerva, need to set internet proxies to download

require(Seurat)
require(R.utils)
require(here)

# CARRARO ALI
# From GEO/Carraro et al: GSE150674_Seurat_Object_ALI.rds.gz
refs <- list(
  "carraro_ALI" = 1,
  "carraro_hBE" = 2
)

fhs <- c(
  here("analysis/supporting_data/cell_references/GSE150674_Seurat_Object_ALI.rds"),
  here("analysis/supporting_data/cell_references/GSE150674_Seurat_Object.rds")
)

https <- c(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE150674&format=file&file=GSE150674%5FSeurat%5FObject%5FALI%2Erds%2Egz",
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE150674&format=file&file=GSE150674%5FSeurat%5FObject%2Erds%2Egz"
)

options(timeout = 3600)
carraro_datasets <- lapply(
  X = refs,
  FUN = function(x) {
    if (!file.exists(fhs[x])) {
      download.file(https[x], paste0(fhs[x], ".gz"), method = "auto")
      gunzip(paste0(fhs[x], ".gz"))
    }
    obj <- readRDS(fhs[x])
    # Subset to only the healthy control samples ("CO")
    Idents(obj) <- "type"
    obj <- subset(obj, idents = "CO")
    DefaultAssay(obj) <- "integrated"
    obj <- UpdateSeuratObject(obj)
    obj[["integrated"]] <- as(object = obj[["integrated"]], Class = "SCTAssay")
    return(obj)
  }
)
options(timeout = 60)

