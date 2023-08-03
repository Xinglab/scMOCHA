# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Aug  3 16:03:51 2023
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

meta_syn12514624 <- data.table::fread(
  input = "/scr1/users/liuc9/mitochondrial/realdata/03-ADKP/SYNAPSE_METADATA_MANIFEST.tsv",
  sep = "\t"
)
# body --------------------------------------------------------------------
meta_syn12514624 |> 
  dplyr::glimpse()




# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------