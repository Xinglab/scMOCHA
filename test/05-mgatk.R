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
  ) ->
  tileplot

ggsave(
  filename = "mgatk_cell_genotype_tileplot.pdf",
  plot = tileplot,
  device = "pdf",
  path = "data/PBMC_10k_v3_10x/result/04-allele-freq/",
  width = 12,
  height = 9
)

depthtable <- readr::read_tsv(file = "/home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/MTbam/mgatk/final/mgatk.depthTable.txt", col_names = c("barcode", "depth"))

depthtable %>% 
  ggplot(aes(x = depth)) +
  geom_density()

depthtable$depth %>% quantile()


coverage <- vroom::vroom(file = "/scr1/users/liuc9/mitochondrial/PBMC_10k_v3_10x/MTbam/mgatk/final/mgatk.coverage.txt.gz", delim = ",", col_names = c("pos", "barcode", "depth"))

coverage %>% 
  tidyr::pivot_wider(
    names_from = pos,
    values_from = depth
  ) ->
  coverage_wider

# seurat ------------------------------------------------------------------



sct_cluster <- readr::read_rds("data/PBMC_10k_v3_10x/rda/pbmc_sct_cluster_annotated.rds.gz")


umap <- sct_cluster@reductions$umap@cell.embeddings
colnames(umap) <- c("UMAP_1", "UMAP_2")

cluster <- sct_cluster@meta.data[, c("seurat_clusters", "sctype")] %>% 
  dplyr::rename(cluster = seurat_clusters, celltype =  sctype)

hetero %>% 
  tidyr::pivot_wider(names_from  = variant, values_from = af) ->
  hetero_w

dplyr::bind_cols(umap, cluster) %>% 
  tibble::rownames_to_column(var = "barcode") %>% 
  dplyr::left_join(hetero_w, by = "barcode") %>% 
  dplyr::mutate(
    dplyr::across(
    dplyr::everything(), 
    ~tidyr::replace_na(.x, 0)
    )) ->
  cell_cluster_af

cell_cluster_af %>% 
  dplyr::left_join(depthtable, by = "barcode") ->
  cell_cluster_af_depth

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
myPalette <- colorRampPalette(c("#969696", "#fa0202"))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 1))

cell_cluster_af %>% 
  # dplyr::filter(variant == "2617A>G") %>% nrow()
  ggplot() +
  geom_point(
    aes(
      x = UMAP_1,
      y = UMAP_2,
      colour = `9966G>A`,
      shape = NULL,
      alpha = NULL
    ),
    size = 0.7
  ) +
  sc


cell_cluster_af_depth %>% 
  ggplot(aes(
    x = depth,
    y = `9966G>A`,
  )) +
  geom_point()


cell_cluster_af_depth %>% 
  dplyr::select(barcode, depth, `9966G>A`) %>% 
  dplyr::inner_join(
    coverage_wider %>% 
      dplyr::select(barcode, `9966`),
    by = "barcode"
  ) %>% 
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(), 
      ~tidyr::replace_na(.x, 0)
    )) %>% 
  ggplot(aes(
    x = `9966`,
    y = `9966G>A`,
  )) +
  geom_point(position = position_jitter())

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(
  file = "data/PBMC_10k_v3_10x/rda/05-mgatk.rda"
)
