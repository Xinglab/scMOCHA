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
  dplyr::mutate(cluster = glue::glue("cluster_{cluster_label}")) |> 
  dplyr::mutate(cluster = forcats::fct_reorder(cluster, cluster_label)) |> 
  dplyr::select(sample_id, cluster) |> 
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
  Seurat::FindClusters(resolution = 0.2) ->
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



sctu_tsne <- sctu |> 
  Seurat::RunTSNE(dims = 1:30, return.model = TRUE) 

FeaturePlot(
  object = sctu_tsne,
  features = c("PTPRC"),
  cols = c("grey", "gold", "#F02415"),
  order = TRUE,
  reduction = "tsne"
) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  ) ->
  p1

FeaturePlot(
  object = sctu_tsne,
  features = c("ISG15"),
  cols = c("grey", "gold", "#F02415"),
  order = TRUE,
  reduction = "tsne"
) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  ) ->
  p2

FeaturePlot(
  object = sctu_tsne,
  features = c("CD83"),
  cols = c("grey", "gold", "#F02415"),
  order = TRUE,
  reduction = "tsne"
) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  ) ->
  p3

FeaturePlot(
  object = sctu_tsne,
  features = c("CD74"),
  cols = c("grey", "gold", "#F02415"),
  order = TRUE,
  reduction = "tsne"
) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  ) ->
  p4


p <- (p1|p2)/(p3|p4) +plot_layout(guides = 'collect')

p
ggsave(
  filename = "NatCommupaper-ref-dotplot.pdf",
  plot = p,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA/03-ADKP/azimuth_motorcortex",
  width = 8,
  height = 7
)

sctu_tsne@assays
DefaultAssay(sctu_tsne) <- "SCT"

FeaturePlot(
  object = sctu_tsne,
  features = c("CD74"),
  cols = c("grey", "gold", "#F02415"),
  order = TRUE,
  reduction = "tsne"
) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )


sctu_tsne$seurat_clusters |> 
  table() ->
  cl

ccl <- factor(glue::glue("{names(cl)} (n={cl})"), glue::glue("{names(cl)} (n={cl})"))

sctu_tsne$seurat_clusters_n <- ccl[sctu_tsne$seurat_clusters]


DimPlot(
  object = sctu_tsne,
  reduction = "tsne",
  cols = paletteer::paletteer_d(
    palette = "ggsci::springfield_simpsons",
    direction = -1
  ),
  group.by = "seurat_clusters_n"
) ->
  pp

ggsave(
  filename = "NatCommupaper-ref-cluster.pdf",
  plot = pp,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA/03-ADKP/azimuth_motorcortex",
  width = 8,
  height = 7
)


Idents(sctu_tsne) |> 
  table() ->
  cl

ccl <- factor(glue::glue("{names(cl)} (n={cl})"), glue::glue("{names(cl)} (n={cl})"))

sctu_tsne$Olha <- ccl[Idents(sctu_tsne)]

DimPlot(
  object = sctu_tsne,
  reduction = "tsne",
  cols = paletteer::paletteer_d(
    palette = "ggsci::springfield_simpsons",
    direction = -1
  ),
  group.by = "Olha"
) ->
  ppp

ggsave(
  filename = "NatCommupaper-ref-cluster-olha.pdf",
  plot = ppp,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA/03-ADKP/azimuth_motorcortex",
  width = 9,
  height = 7
)



# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------