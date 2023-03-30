# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Mar 30 18:05:33 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(Seurat)
library(Azimuth)


# args --------------------------------------------------------------------

args <- commandArgs(TRUE)

h5file <- args[1]
refname <- args[2]
# h5file <- "/scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/flu2/Flu2/outs/filtered_feature_bc_matrix.h5"
# refname <- "pbmcref"

# src ---------------------------------------------------------------------

pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv")


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------


# body --------------------------------------------------------------------

sc <- Azimuth::RunAzimuth(
  query = h5file,
  reference = refname
)


# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------