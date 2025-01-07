#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Fri May 24 14:54:13 2024
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
library(prismatic)
library(paletteer)
library(data.table)
#library(rlang)
library(GetoptLong)
library(logger)

# args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.
verbose <- FALSE
spec <- "
Usage: Rscript foorbar.R [options]

Options:
<verbose!> Print messages
"

GetoptLong.options(help_style = "two-column")
GetoptLong(spec, template_control = list(opt_width = 21))

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------
log_threshold(TRACE)
log_layout(layout_glue_colors)

# log_info('Starting the script...')
# log_debug('This is the second log line')
# log_trace('Note that the 2nd line is being placed right after the 1st one.')
# log_success('Doing pretty well so far!')
# log_warn('But beware, as some errors might come :/')
# log_error('This is a problem')
# log_debug('Note that getting an error is usually bad')
# log_error('This is another problem')
# log_fatal('The last problem')

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------

d <- readr::read_csv("/scr1/users/liuc9/mitochondrial/realdata/05-Liming/disco/DISCO_Marker_blood.csv")

# body --------------------------------------------------------------------

d |> 
  dplyr::select(
    disease,
    sampleId,
    sampleType,
    platform,
    projectId,
    age
  ) ->
  dd

dd |> 
  dplyr::mutate(
    disease = ifelse(is.na(disease), "Normal", disease)
  ) |> 
  dplyr::count(sampleType, disease) |> 
  dplyr::mutate(disease_r = n / sum(n)) |> 
  dplyr::group_by(sampleType) |> 
  dplyr::mutate(sampleType_r = sum(disease_r)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(
    sampleType = "{sampleType} {round(sampleType_r * 100, 2)}%" |> glue::glue(),
    disease = "{disease} {round(disease_r * 100, 2)} %" |> glue::glue()
  ) |> 
  dplyr::select(1, 2, 3) |> 
  plotme::count_to_sunburst()


dd |> 
  dplyr::count(projectId) |> 
  dplyr::arrange(-n) |> 
  dplyr::filter(grepl("GSE", projectId)) |> 
  dplyr::filter(n >5) ->
  dd_n_samples

dd |> 
  dplyr::filter(projectId %in% dd_n_samples$projectId) |> 
  dplyr::select(disease, projectId) |> 
  dplyr::count(projectId, disease)


# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------