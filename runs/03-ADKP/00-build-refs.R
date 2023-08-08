#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Aug  8 14:14:19 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
#library(rlang)
library(Matrix)
library(Seurat)
library(presto)
library(Azimuth)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
annotations <- readr::read_csv(
  file = "/scr1/users/liuc9/mitochondrial/realdata/03-ADKP/forrefs/annofile.csv"
) |> 
  dplyr::select(sample_id, cluster_label) |> 
  tibble::deframe()
# countmatrix <- vroom::vroom(
#   file = "/scr1/users/liuc9/mitochondrial/realdata/03-ADKP/forrefs/countmatrix.csv",
#   delim = ","
# ) |> 
#   dplyr::rename(genename = "...1") 

countmatrix <- data.table::fread(
  input = "/scr1/users/liuc9/mitochondrial/realdata/03-ADKP/forrefs/countmatrix.csv",
  sep = ","
) |> 
  tibble::column_to_rownames(var = "V1")

countmatrix[1:4, 1:4]


# body --------------------------------------------------------------------

sc <- Seurat::CreateSeuratObject(
  counts = countmatrix,
  project = "syn21438358"
)

sct <- Seurat::SCTransform(object = sc)

sct |>
  Seurat::RunPCA() |>
  Seurat::RunUMAP(dims = 1:30, return.model = TRUE) |>
  Seurat::FindNeighbors(dims = 1:30) |>
  Seurat::FindClusters(resolution = c(0.8, 2)) ->
  sctu

annotations
Idents(sctu) |>
  names() ->
  cells

cells_anno <- annotations[cells] |> tidyr::replace_na("remove")
names(cells_anno) <- cells
cells_anno <- factor(cells_anno)
Idents(object = sctu) <- cells_anno


sctu <- RenameCells(
  object = sctu,
  new.names = unname(obj = sapply(
    X = Seurat::Cells(x = sctu),
    FUN = function(.s) {
      paste0("syn21438358", .s)
    }
    
  ))
)

ref <- sctu

if ("remove" %in% levels(x = ref)) {
  ref <- subset(x = ref, idents = "remove", invert = TRUE)
  ref <- RunPCA(object = ref, verbose = FALSE)
}
ref$annotation.l1 <- Idents(object = ref)
ref <- RunUMAP(object = ref, dims = 1:30, return.model = TRUE)
full.ref <- ref
colormap <- list(annotation.l1 = CreateColorMap(object = ref, seed = 2))
colormap[["annotation.l1"]] <- colormap[["annotation.l1"]][sort(x = names(x = colormap[["annotation.l1"]]))]

ref <- AzimuthReference(
  object = ref,
  refUMAP = "umap",
  refDR = "pca",
  refAssay = "SCT",
  metadata = c("annotation.l1"),
  dims = 1:50,
  k.param = 31,
  colormap = colormap,
  reference.version = "1.0.0"
)

ref.dir <- "/home/liuc9/github/scMOCHA/03-ADKP/forrefs/azimuth_syn21438358"
dir.create(
  path = ref.dir,
  showWarnings = F,
  recursive = T
)

SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]], file = file.path(ref.dir, "idx.annoy"))
saveRDS(object = ref, file = file.path(ref.dir, "ref.Rds"))
saveRDS(object = full.ref, file = file.path(ref.dir, "fullref.Rds"))

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------