# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Fri Feb 24 01:59:16 2023
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

rootdir <- "/home/liuc9/scratch/mitochondrial/testdata/2_smart_seq2"

sratable <- readr::read_csv(
  file = file.path(
    rootdir,
    "SraRunTable.txt"
  )
)
# body --------------------------------------------------------------------

sratable %>%
  dplyr::mutate(prefetch = "prefetch --max-size 50G {Run} --output-directory {rootdir} &" %>% glue::glue()) %>%
  dplyr::select(Bytes, prefetch) %>%
  dplyr::arrange(Bytes) %>%
  dplyr::select(prefetch) %>% 
  readr::write_tsv(file = "{rootdir}/sratable_prefetch.sh" %>% glue::glue(), col_names = F)

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------