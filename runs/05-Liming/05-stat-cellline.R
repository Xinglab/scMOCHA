#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed Dec 13 15:24:38 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
library(data.table)
#library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
datadir <- "/home/liuc9/github/scMOCHA/05-Liming/cellline/torun"

# body --------------------------------------------------------------------

tibble::tibble(
  torun = list.dirs(
    path = datadir,
    full.names = T,
    recursive = F 
  )
) |> 
  dplyr::mutate(
    log = purrr::map_chr(
      .x = torun,
      .f = \(.x) {
        .log <- "{basename(.x)}.log" |> glue::glue()
        file.path(
          .x,
          .log
        )
      }
    )
  ) |> 
  dplyr::mutate(
    projectname = basename(torun)
  ) ->
  logfile

logfile |> 
  dplyr::mutate(
    outputdir = purrr::map_chr(
      .x = log,
      .f = \(.log) {
        # .log <- logfile$log[[1]]
        
        .l <- readr::read_lines(
          file = .log
        ) 
        
        .l |> 
          stringr::str_detect(
            "scMOCHA.output_dir_tar_gz"
          ) |> 
          which() ->
          .arr
        
        if(length(.arr) == 0) {return(NA_character_)}
        .l[[.arr[[1]]]] |> 
          gsub("\"| |,|.tar.gz", "", x = _) |> 
          strsplit(
            split = ":"
          ) ->
          .s
        .s[[1]][[2]]
      }
    )
  ) |> 
  dplyr::select(
    projectname, outputdir
  ) |> 
  dplyr::filter(!is.na(outputdir)) ->
  outdir

outdir |> 
  dplyr::mutate(
    cluster = purrr::map(
      .x = outputdir,
      .f = \(.x) {
        # .x <- outdir$outputdir[[1]]
        qc_cell_stats <- readxl::read_xlsx(
          path = file.path(
            .x, 
            "qc_cell_stats.xlsx"
          )
        )
        read_depth <- data.table::fread(
          input = file.path(
            .x,
            "possorted_genome_bam.MT.depth"
          )
        )
        celltype_ratio <- readr::read_tsv(
          file.path(
            .x,
            "celltype_ratio.tsv"
          )
        )
        cell_variant_annotation <- readr::read_tsv(
          file.path(
            .x,
            "cell_variant_annotation.tsv"
          )
        )
        
        cell_heteroplasmic_df_raw <- readr::read_tsv(
          file.path(
            .x,
            "cell.cell_heteroplasmic_df_raw.tsv.gz"
          )
        )
        cell_coverage <- readr::read_tsv(
          file.path(
            .x,
            "cell.coverage.txt.gz"
          )
        )
        
        tibble::tibble(
          qc_cell_stats = list(qc_cell_stats),
          read_depth = list(read_depth),
          celltype_ratio = list(celltype_ratio),
          cell_variant_annotation = list(cell_variant_annotation),
          cell_heteroplasmic_df_raw = list(cell_heteroplasmic_df_raw),
          cell_coverage = list(cell_coverage)
        )
      }
    )
  ) |> 
  tidyr::unnest(cols = a) ->
  alldataloaded

alldataloaded


# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(
  "/home/liuc9/github/scMOCHA/05-Liming/cellline/torun/05-stat-cellline.rda"
)