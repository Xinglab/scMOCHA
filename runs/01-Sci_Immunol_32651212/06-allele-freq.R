# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu May  4 15:46:02 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
metadata_anno_depth <- 
  readr::read_rds(
    file = "/home/liuc9/github/scRNAseq-MitoVariant/01-Sci_Immunol_32651212/outputs/metadata_anno_depth.rds"
  )


# body --------------------------------------------------------------------

metadata_anno_depth |> dplyr::glimpse()

metadata_anno_depth |> 
  dplyr::mutate(
    variant = purrr::map2(
      .x = anno,
      .y = tardir,
      .f = function(.x, .y) {
        if(is.na(.y)) {return(NULL)}
        .x |> 
          dplyr::mutate(
            variant = glue::glue("{tpos}{tnt}>{qnt}")
          ) |> 
          dplyr::pull(variant)
      }
    )
  ) ->
  metadata_anno_depth_variant


metadata_anno_depth_variant |> 
  dplyr::select(srrid, source_name, variant) |> 
  dplyr::filter(!purrr::map_lgl(.x = variant, .f = is.null)) ->
  for_variant


for_variant |> 
  dplyr::filter(source_name == "Normal_PBMC") |> 
  dplyr::select(srrid, variant) |> 
  tibble::deframe() ->
  d

ggvenn::ggvenn(
  data = d
)
# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------