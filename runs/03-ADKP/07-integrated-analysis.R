#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Fri Aug 18 18:11:28 2023
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

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------

outfiles <- readr::read_tsv(
  file = "/home/liuc9/github/scMOCHA/03-ADKP/output/outfiles.tsv"
)


 <- readr::read_tsv(
  file = "/home/liuc9/github/scMOCHA/03-ADKP/SYNAPSE_METADATA_MANIFEST.tsv"
) |>
  dplyr::rename(
    srrid = Run
  )

# body --------------------------------------------------------------------



# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------