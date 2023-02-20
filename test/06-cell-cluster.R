# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Mon Feb 20 16:21:10 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
sct_cluster <- readr::read_rds(
  file = "data/PBMC_10k_v3_10x/rda/pbmc_sct_cluster.rds.gz"
)


sct_cluster@meta.data %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(
    var = "cellbarcode"
  ) %>% 
  dplyr::select(
    cellbarcode,
    seurat_clusters
  ) %>% 
  tibble::as_tibble() %>% 
  dplyr::mutate(tag = "CJ") %>% 
  dplyr::select(1, 3, 2) ->
  cellbarcode

# body --------------------------------------------------------------------

cellbarcode %>% 
  readr::write_tsv(
    file = "/home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/barcode_cluster.tsv",
    col_names = F
  )

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------