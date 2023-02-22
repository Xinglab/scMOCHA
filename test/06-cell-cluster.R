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
  dplyr::mutate(seurat_clusters = as.numeric(seurat_clusters)) %>% 
  dplyr::mutate(
    cluster = purrr::map_chr(
      .x = seurat_clusters,
      .f = function(.x) {
        # if (.x < 9) {
        #   glue::glue("GATTACAAcluster{.x}-1")
        # } else(
        #   glue::glue("GATTACAcluster{.x}-1")
        # )
        glue::glue("cluster{.x}-1")
      }
    )
  ) %>% 
  dplyr::select(1, 3, 4) ->
  cellbarcode



# body --------------------------------------------------------------------

cellbarcode %>% 
  readr::write_tsv(
    file = "/home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/barcode_cluster.tsv",
    col_names = F
  )

cellbarcode %>% 
  dplyr::select(cluster) %>% 
  dplyr::arrange(cluster) %>% 
  dplyr::distinct() %>% 
  readr::write_tsv(
    file = "/home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/cluster.tsv",
    col_names = F
  )

cellbarcode %>% 
  dplyr::mutate(cluster = "Bulk") %>% 
  readr::write_tsv(
    file = "/home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/barcode_bulk.tsv",
    col_names = F
  )
  

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------