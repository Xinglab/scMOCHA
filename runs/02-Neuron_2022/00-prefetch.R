# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Mon May  1 22:30:51 2023
# @DESCRIPTION: 


# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------


# body --------------------------------------------------------------------
sratable <- readr::read_csv(
  file = "/scr1/users/liuc9/mitochondrial/realdata/02-Neuron_2022/SraRunTable.txt"
)


sratable %>%
  dplyr::mutate(prefetch = "prefetch --max-size 50G {Run} --output-directory /scr1/users/liuc9/mitochondrial/realdata/02-Neuron_2022/ &" %>% glue::glue()) %>%
  dplyr::select(Bytes, prefetch) %>%
  dplyr::arrange(Bytes) %>%
  dplyr::select(prefetch) %>% 
  readr::write_tsv(file = "/scr1/users/liuc9/mitochondrial/realdata/02-Neuron_2022/sratable_prefetch.sh", col_names = F)

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------