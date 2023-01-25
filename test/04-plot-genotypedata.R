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
  file = "data/PBMC_10k_v3_10x/test/cellsnp/cellSNP.cells.vcf.gz",
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
  dplyr::mutate(gt = furrr::future_map(
    .x = genotype,
    .f = function(.x) {
      a <- strsplit(.x, split=":")[[1]]
      # print(a)
      gt <- a[[1]]
      ad <- a[[2]]
      dp <- a[[3]]
      
      tibble::tibble(
        gt = gt,
        ad = ad,
        dp = dp
      )
    }
  )) %>% 
  tidyr::unnest(cols = gt) ->
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
  dplyr::mutate(
    ad = ifelse(ad == ".", 0, as.integer(ad))
  ) %>% 
  dplyr::mutate(
    dp = ifelse(dp == ".", 0, as.integer(dp))
  ) %>% 
  # dplyr::select(variant, barcode, a) %>% 
  dplyr::mutate(a = factor(a)) %>% 
  dplyr::mutate(af = ad / dp) ->
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
  p;p

ggsave(
  filename = "cell_genotype.pdf",
  plot = p,
  device = "pdf",
  path = "data/PBMC_10k_v3_10x/result/02-cluster",
  width = 12,
  height = 9
)


# af ----------------------------------------------------------------------

cell_variant_gt %>% 
  dplyr::filter(ad > 1) %>% 
  dplyr::mutate(af = ad / dp) %>% 
  dplyr::filter(!is.nan(af)) ->
  cell_variant_gt_af

cell_variant_gt_af %>% 
  dplyr::group_by(variant) %>% 
  dplyr::summarise(s_af = sum(af)) %>% 
  dplyr::arrange(s_af) ->
  cell_variant_gt_af_variant_rank

cell_variant_gt_af %>% 
  dplyr::group_by(barcode) %>% 
  dplyr::summarise(s_af = sum(af)) %>% 
  dplyr::arrange(-s_af) ->
  cell_variant_gt_af_cell_rank
CPCOLS <- c("#191970", "#F8F8FF", "#FF4040")
cell_variant_gt_af %>% 
  dplyr::mutate(variant = factor(variant, levels = cell_variant_gt_af_variant_rank$variant)) %>% 
  dplyr::mutate(barcode = factor(barcode, levels = cell_variant_gt_af_cell_rank$barcode)) %>% 
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
  )  ->
  tileplot;tileplot
ggsave(
  filename = "cell_genotype_tileplot.pdf",
  plot = tileplot,
  device = "pdf",
  path = "data/PBMC_10k_v3_10x/result/02-cluster",
  width = 12,
  height = 9
)

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------

save.image(file = "data/PBMC_10k_v3_10x/rda/04-plot-genotypedata.rda")
