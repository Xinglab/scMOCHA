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


# load data --------------------------------------------------------------

args <- commandArgs(TRUE)

hetero_file <- args[1]
coverage_file <- args[2]
cluster_umap_file <- args[3]

# hetero_file <- "/home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/MTbam/mgatk/final/mgatk.cell_heteroplasmic_df.tsv.gz"

hetero <- vroom::vroom(file = hetero_file) %>%
  dplyr::rename(barcode = `...1`) %>%
  tidyr::pivot_longer(
    cols = -barcode,
    names_to = "variant",
    values_to = "af"
  )


coverage <- vroom::vroom(file = coverage_file, delim = ",", col_names = c("pos", "barcode", "depth"))

coverage %>%
  dplyr::mutate(depth = log2(depth + 1)) ->
  coverage_log2

coverage_log2 %>%
  tidyr::pivot_wider(
    names_from = pos,
    values_from = depth
  ) ->
  coverage_wider

cluster_umap <- vroom::vroom(file = cluster_umap_file)


# body --------------------------------------------------------------------

hetero %>%
  dplyr::group_by(variant) %>%
  dplyr::summarise(s_af = sum(af, na.rm = T)) %>%
  dplyr::arrange(s_af) ->
  hetero_variant_rank

hetero %>%
  dplyr::group_by(barcode) %>%
  dplyr::summarise(s_af = sum(af)) %>%
  dplyr::arrange(s_af) ->
  hetero_barcode_rank


# seurat ------------------------------------------------------------------


hetero %>%
  tidyr::pivot_wider(names_from  = variant, values_from = af) ->
  hetero_w

cluster_umap %>%
  tibble::rownames_to_column(var = "barcode") %>%
  dplyr::left_join(hetero_w, by = "barcode") ->
  cell_cluster_af




# Heatmap -----------------------------------------------------------------



cell_cluster_af %>%
  dplyr::select(barcode, dplyr::contains(">")) %>%
  tidyr::pivot_longer(
    cols = -barcode,
    names_to = "variant",
    values_to = "af"
  ) %>%
  dplyr::group_by(barcode) %>%
  dplyr::summarise(s_af = sum(af, na.rm = T)) %>%
  dplyr::left_join(
    cell_cluster_af %>%
      dplyr::select(barcode, cluster) ,
    by = "barcode"
  ) %>%
  dplyr::arrange(cluster, -s_af) ->
  cell_cluster_af_col_rank

cell_cluster_af %>%
  dplyr::select(barcode, dplyr::contains(">")) %>%
  tidyr::pivot_longer(
    cols = -barcode,
    names_to = "variant",
    values_to = "af"
  ) %>%
  dplyr::mutate(pos = gsub(pattern = "([[:digit:]]*).*", "\\1", variant) %>% as.numeric()) ->
  cell_cluster_af_pos


coverage_log2 %>%
  dplyr::filter(barcode %in% cell_cluster_af_pos$barcode) %>%
  dplyr::filter(pos %in% cell_cluster_af_pos$pos) ->
  coverage_log2_pos

cell_cluster_af_pos %>%
  dplyr::left_join(
    coverage_log2_pos,
    by = c("barcode", "pos")
  ) %>%
  tidyr::replace_na(
    replace = list(
      af = 0
    )
  ) %>%
  dplyr::mutate(af = ifelse(is.na(depth), NA, af)) %>%
  dplyr::arrange(pos) %>%
  dplyr::select(barcode, variant, af) %>%
  tidyr::pivot_wider(
    names_from = "variant",
    values_from = af
  ) %>%
  dplyr::slice(match(cell_cluster_af_col_rank$barcode, barcode)) %>%
  tibble::column_to_rownames(var = "barcode") %>%
  as.matrix() %>%
  t() ->
  cell_cluster_af_mtx

library(ComplexHeatmap)
pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv")


cell_cluster_af %>%
  dplyr::select(barcode, cluster, celltype ) %>%
  dplyr::slice(match(colnames(cell_cluster_af_mtx), barcode)) %>%
  tibble::column_to_rownames(var = "barcode") %>%
  dplyr::select(Cluster = cluster) ->
  cell_cluster_af_cluster

col_clusters <- levels(cell_cluster_af_cluster$Cluster) %>% as.numeric()
col_colors <- pcc$color[1:length(levels(cell_cluster_af_cluster$Cluster))]

names(col_colors) <- col_clusters

chm_top <- ComplexHeatmap::HeatmapAnnotation(
  df = cell_cluster_af_cluster,
  gap = unit(c(2,2), "mm"),
  col = list(Cluster = col_colors),
  which = "column"
)


ComplexHeatmap::Heatmap(
  matrix = cell_cluster_af_mtx,
  col = circlize::colorRamp2(
    breaks = c(0, 1),
    colors = c("white", "red"),
    space = "RGB"
  ),
  name = "Allele Freq",
  na_col = "grey",
  color_space = "LAB",
  rect_gp = gpar(col = NA),
  border = NA,
  cell_fun = NULL,
  layer_fun = NULL,
  jitter = FALSE,
  # row
  cluster_rows = F,
  cluster_row_slices = T,
  clustering_distance_rows = "pearson",
  clustering_method_rows = "ward.D",
  # column
  cluster_columns = FALSE,
  cluster_column_slices = T,
  # clustering_distance_columns = "pearson",
  # clustering_method_columns = "ward.D",
  show_column_names = FALSE,

  top_annotation = chm_top
) ->
  chm;chm


{
  pdf(
    file = "mgatk_cell_genotype_heatmap.pdf",
    width = 12,
    height = 4
  )
  ComplexHeatmap::draw(object = chm)

  dev.off()
}

cell_cluster_af_pos %>%
  dplyr::left_join(
    coverage_log2_pos,
    by = c("barcode", "pos")
  ) %>%
  dplyr::select(barcode, pos, depth) %>%
  dplyr::arrange(pos) %>%
  tidyr::pivot_wider(
    names_from = pos,
    values_from = depth
  ) %>%
  dplyr::slice(match(cell_cluster_af_col_rank$barcode, barcode)) %>%
  tibble::column_to_rownames(var = "barcode") %>%
  as.matrix() %>%
  t() ->
  depth_mtx

ComplexHeatmap::Heatmap(
  matrix = depth_mtx,
  col = circlize::colorRamp2(
    breaks = c(1, 4),
    colors = c("white", "red"),
    space = "RGB"
  ),
  name = "log2(Depth+1)",
  na_col = "grey",
  color_space = "LAB",
  rect_gp = gpar(col = NA),
  border = NA,
  cell_fun = NULL,
  layer_fun = NULL,
  jitter = FALSE,
  # row
  cluster_rows = F,
  cluster_row_slices = T,
  clustering_distance_rows = "pearson",
  clustering_method_rows = "ward.D",
  # column
  cluster_columns = FALSE,
  cluster_column_slices = T,
  # clustering_distance_columns = "pearson",
  # clustering_method_columns = "ward.D",
  show_column_names = FALSE,

  top_annotation = chm_top
) ->
  chm_depth;chm_depth

{
  pdf(
    file = "mgatk_cell_depth_heatmap.pdf",
    width = 12,
    height = 4
  )
  ComplexHeatmap::draw(object = chm_depth)

  dev.off()
}


# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
# save.image(file = "data/PBMC_10k_v3_10x/rda/05-mgatk.rda")

# load(file = "data/PBMC_10k_v3_10x/rda/05-mgatk.rda")
