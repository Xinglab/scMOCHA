# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Dec  1 11:34:25 2022
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

# load vcf file
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)
vcf <- readr::read_tsv(
  file = "data/test/cellsnp/cellSNP.cells.vcf.gz",
  comment = "##"
)



# body --------------------------------------------------------------------

vcf %>% 
  dplyr::select(1:9) %>% 
  dplyr::mutate(a = glue::glue("{REF}{POS}{ALT}")) ->
  vcf_variant
vcf %>% 
  dplyr::select(-c(1:9)) %>% 
  tibble::add_column(variant = vcf_variant$a) %>% 
  tidyr::gather(key = barcode, value = genotype, -variant) %>% 
  dplyr::mutate(gt = purrr::map_chr(
    .x = genotype,
    .f = function(.x) {
      strsplit(.x, split=":")[[1]][1]
    }
  )) ->
  cell_variant

cell_variant %>% 
  dplyr::mutate(
    a = dplyr::case_when(
      gt == "." ~ 0,
      gt == "0/0" ~ 1,
      gt == "1/0" ~ 2,
      gt == "0/1" ~ 2,
      gt == "1/1" ~ 3
    )
  ) %>% 
  dplyr::select(variant, barcode, a) %>% 
  dplyr::mutate(a = factor(a)) ->
  cell_variant_gt

cell_variant_gt %>% 
  dplyr::group_by(variant) %>% 
  dplyr::summarise(s = sum(as.numeric(a))) %>% 
  dplyr::arrange(s) ->
  variant_rank

variant_rank %>% 
  dplyr::arrange(-s) %>% 
  head(10) %>% 
  dplyr::pull(variant)

cell_variant_gt %>% 
  dplyr::mutate(variant = factor(variant, levels = variant_rank$variant)) %>% 
  ggplot(aes(x = barcode, y = variant)) +
  geom_tile(aes(fill = a)) +
  scale_fill_manual(
    name = "Genotype",
    labels = c("", "0/0", "1/0", "1/1"),
    values = c("white", ggsci::pal_aaas()(3))
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    x = "Cell",
    y = "MT variant"
  ) ->
  p

ggsave(
  filename = "cell_genotype.pdf",
  plot = p,
  device = "pdf",
  path = "data/result/02-cluster",
  width = 12,
  height = 9
)


# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------

save.image(file = "data/rda/04-plot-genotypedata.rda")
