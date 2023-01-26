# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Jan 26 00:13:25 2023
# @DESCRIPTION: 05-mgatk.R

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
hetero_file <- "/home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/MTbam/mgatk/final/mgatk.cell_heteroplasmic_df.tsv.gz"

hetero <- vroom::vroom(file = hetero_file) %>% 
  dplyr::rename(barcode = `...1`) %>% 
  tidyr::pivot_longer(
    cols = -barcode,
    names_to = "variant",
    values_to = "af"
  ) %>% 
  tidyr::replace_na(list(af = 0))


# body --------------------------------------------------------------------

hetero %>% 
  dplyr::group_by(variant) %>% 
  dplyr::summarise(s_af = sum(af)) %>% 
  dplyr::arrange(s_af) ->
  hetero_variant_rank

hetero %>% 
  dplyr::group_by(barcode) %>% 
  dplyr::summarise(s_af = sum(af)) %>% 
  dplyr::arrange(s_af) ->
  hetero_barcode_rank

hetero %>% 
  dplyr::mutate(variant = factor(variant, levels = hetero_variant_rank$variant)) %>% 
  dplyr::mutate(barcode = factor(barcode, levels = hetero_barcode_rank$barcode)) %>% 
  dplyr::rename(`Allele Freq` = af) %>% 
  ggplot(aes(
    x = barcode,
    y = variant
  )) +
  geom_tile(aes(fill = `Allele Freq`)) +
  # scale_fill_gradient2(
  #   low = "white",
  #   mid = CPCOLS[3],
  #   high = "#67000e"
  # ) +
  # scale_fill_brewer(name = "Alelle Freq") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    x = "Cell",
    y = "MT variant"
  )


# seurat ------------------------------------------------------------------



sct_cluster <- readr::read_rds()

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------