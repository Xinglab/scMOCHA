#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Aug 10 15:34:40 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
#library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------


# body --------------------------------------------------------------------
SeuratData::AvailableData() |> 
  dplyr::filter(grepl("Azimuth Reference", x = Summary)) |> 
  dplyr::select(system, Dataset, species, ncells, tech) |> 
  data.table::as.data.table() |> 
  dplyr::mutate(
    system = purrr::map_chr(
      .x = system, 
      .f = stringr::str_to_title
    )
  ) |> 
  dplyr::arrange(-ncells) ->
  azd

# SeuratData::LoadData(
#   ds = azd$Dataset[[1]],
#   type = "azimuth"
# )
# 
# .a <- .Last.value
# .a$plot@meta.data |> 
#   dplyr::select(-c(orig.ident, nCount_RNA, nFeature_RNA)) |> 
#   data.table::as.data.table() |> 
#   dplyr::distinct() 
{
  future::plan(future::multisession, workers = 10)
  azd |> 
    dplyr::mutate(
      a = furrr::future_map(
        .x = Dataset,
        .f = \(.ds) {
          .d <- SeuratData::LoadData(
            ds = .ds,
            type = "azimuth"
          )
          .d$plot@meta.data |>
            dplyr::select(-c(orig.ident, nCount_RNA, nFeature_RNA)) |>
            data.table::as.data.table() |>
            dplyr::distinct()
        }
      )
    ) ->
    azd_celltype
  future::plan(future::sequential)
}



# footer ------------------------------------------------------------------


# save image --------------------------------------------------------------