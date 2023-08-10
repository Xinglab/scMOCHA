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


azd_celltype |> 
  dplyr::mutate(
    b = furrr::future_map_chr(
      .x = a,
      .f = \(.a) {
        as.list(.a) |> 
          purrr::map(.f = levels) |> 
          purrr::map(.f = length) |> 
          tibble::enframe() |> 
          tidyr::unnest(cols = value) |> 
          dplyr::mutate(nv = glue::glue("{name} (n={value})")) ->
          .aa
        
        paste0(.aa$nv, collapse = "; ")
      }
    )
  ) ->
  azd_celltype_n


azd_celltype_n |> 
  dplyr::select(-a) |> 
  dplyr::rename(
    Tissue = system,
    Species = species,
    Ncells = ncells,
    Tech = tech,
    Celltypes = b
  ) |> 
  writexl::write_xlsx(
    path = "/home/liuc9/github/scMOCHA/03-ADKP/azimuth_celltype.xlsx"
  )


# footer ------------------------------------------------------------------


# save image --------------------------------------------------------------