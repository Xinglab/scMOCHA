# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Feb 21 14:22:19 2023
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
hetero <- readr::read_tsv(
  file = "/home/liuc9/tmp/mgatk2/final/mgatk.cell_heteroplasmic_df.tsv.gz"
) %>% 
  dplyr::rename(
    barcode = `...1`
  )

variant <- readr::read_tsv(
  file = "/home/liuc9/tmp/mgatk2/final/mgatk.variant_stats.tsv.gz"
)

variant %>% 
  dplyr::filter(
    n_cells_conf_detected >=3
  ) %>% 
  dplyr::pull(variant) ->
  multi_cell_variants


hetero

# body --------------------------------------------------------------------



# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------