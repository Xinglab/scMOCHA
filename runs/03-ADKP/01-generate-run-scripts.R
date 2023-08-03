# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Aug  3 17:07:54 2023
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

datadir <- "/home/liuc9/github/scMOCHA/03-ADKP"
metadata <- data.table::fread(
  input = file.path(
    datadir,
    "fastq",
    "SYNAPSE_METADATA_MANIFEST.tsv"
  ),
  sep = "\t"
)

metadata |> 
  dplyr::glimpse()

# body --------------------------------------------------------------------

metadata |> 
  dplyr::select(name, specimenID, individualID) |> 
  dplyr::arrange(specimenID, name)

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------