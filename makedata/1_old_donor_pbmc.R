# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Fri Feb 24 01:52:26 2023
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

sratable <- readr::read_csv(
  file = "/scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/SraRunTable.txt"
)

# body --------------------------------------------------------------------
sratable %>%
  dplyr::mutate(prefetch = "prefetch --max-size 50G {Run} --output-directory /scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/ &" %>% glue::glue()) %>%
  dplyr::select(Bytes, prefetch) %>%
  dplyr::arrange(Bytes) %>%
  dplyr::select(prefetch) %>% 
  readr::write_tsv(file = "/scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/sratable_prefetch.sh", col_names = F)


# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------