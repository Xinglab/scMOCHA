# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Feb 21 17:19:24 2023
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

hetero_file <- "/scr1/users/liuc9/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/mgatk_cluster/final/mgatk_cluster.cell_heteroplasmic_df.tsv.gz"
hetero_file_bulk <- "/scr1/users/liuc9/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/mgatk_bulk/final/mgatk_bulk.cell_heteroplasmic_df.tsv.gz"

hetero <- vroom::vroom(file = hetero_file) %>% 
  dplyr::rename(barcode = `...1`) %>% 
  tidyr::pivot_longer(
    cols = -barcode,
    names_to = "variant",
    values_to = "af"
  ) %>% 
  dplyr::mutate(
    barcode = gsub(
      pattern = "cluster|-1",
      replacement = "",
      x = barcode
    )
  ) %>% 
  dplyr::mutate(barcode = as.integer(barcode) -1) 
  # dplyr::mutate(pos = gsub(pattern = "([[:digit:]]*).*", "\\1", variant) %>% as.numeric()) 

hetero_bulk <- vroom::vroom(file = hetero_file_bulk) %>% 
  dplyr::rename(barcode = `...1`) %>% 
  dplyr::mutate(barcode = 100) %>% 
  tidyr::pivot_longer(
    cols = -barcode,
    names_to = "variant",
    values_to = "af"
  ) 

hetero_cb <- dplyr::bind_rows(
  hetero,
  hetero_bulk
) %>% 
  tidyr::pivot_wider(
    names_from = variant,
    values_from = af
  ) %>% 
  tidyr::pivot_longer(
    -barcode,
    names_to = "variant",
    values_to = "af"
  ) %>% 
  dplyr::mutate(pos = gsub(pattern = "([[:digit:]]*).*", "\\1", variant) %>% as.numeric())
  


# depthtable <- readr::read_tsv(file = "/scr1/users/liuc9/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/mgatk_cluster/final/mgatk_cluster.depthTable.txt", col_names = c("barcode", "depth")) %>% 
#   dplyr::mutate(
#     barcode = gsub(
#       pattern = "cluster|-1",
#       replacement = "",
#       x = barcode
#     )
#   ) %>% 
#   dplyr::mutate(barcode = as.integer(barcode) -1)

coverage <- vroom::vroom(file = "/scr1/users/liuc9/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/mgatk_cluster/final/mgatk_cluster.coverage.txt.gz", delim = ",", col_names = c("pos", "barcode", "depth")) %>% 
  dplyr::mutate(
    barcode = gsub(
      pattern = "cluster|-1",
      replacement = "",
      x = barcode
    )
  ) %>% 
  dplyr::mutate(barcode = as.integer(barcode) -1)

coverage_bulk <- vroom::vroom(file = "/scr1/users/liuc9/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/mgatk_bulk/final/mgatk_bulk.coverage.txt.gz", delim = ",", col_names = c("pos", "barcode", "depth")) %>% 
  dplyr::mutate(barcode = 100)

coverage_cb <- dplyr::bind_rows(
  coverage,
  coverage_bulk
) %>% 
  dplyr::filter(pos %in% hetero_cb$pos)


# body --------------------------------------------------------------------



hetero_cb %>% 
  dplyr::inner_join(
    coverage_cb,
    by = c("barcode", "pos")
  ) ->
  hetero_af_pos

hetero_af_pos %>% 
  dplyr::group_by(
    barcode
  ) %>% 
  dplyr::summarise(
    s_af = sum(af, na.rm = T)
  ) %>% 
  dplyr::arrange(barcode) -> 
  hetero_af_pos_rank

hetero_af_pos %>% 
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
  dplyr::slice(match(hetero_af_pos_rank$barcode, barcode)) %>% 
  tibble::column_to_rownames(var = "barcode") %>% 
  as.matrix() %>% 
  t() ->
  cluster_af_mtx



library(ComplexHeatmap)
pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv")

hetero_af_pos_rank %>% 
  dplyr::select(-s_af) %>% 
  dplyr::mutate(cluster = barcode) %>% 
  dplyr::mutate(cluster = factor(cluster)) %>% 
  tibble::column_to_rownames(var = "barcode") %>% 
  dplyr::select(`Cell cluster` = cluster) ->
  cluster_cluster

col_clusters <- levels(cluster_cluster$`Cell cluster`) %>% as.numeric()
col_colors <- pcc$color[1:length(levels(cluster_cluster$`Cell cluster`))]

names(col_colors) <- col_clusters

chm_top <- ComplexHeatmap::HeatmapAnnotation(
  df = cluster_cluster,
  gap = unit(c(2,2), "mm"),
  col = list(`Cell cluster` = col_colors),
  which = "column"
)
c("#FFAEB9", "#FFFFFF", "#FFFFFF")

ComplexHeatmap::Heatmap(
  matrix = cluster_af_mtx,
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
    file = "data/PBMC_10k_v3_10x/result/04-allele-freq/mgatk_cluster_genotype_heatmap-bulk.pdf",
    width = 7, 
    height = 7
  )
  ComplexHeatmap::draw(object = chm)
  
  dev.off()
}

hetero_af_pos %>% 
  dplyr::select(barcode, variant, pos, depth)  %>% 
  dplyr::mutate(depth = log2(depth +1)) ->
  hetero_af_pos_d

hetero_af_pos_d %>% 
  dplyr::arrange(pos) %>% 
  dplyr::select(-pos) %>% 
  tidyr::pivot_wider(
    names_from = variant,
    values_from = depth
  ) %>% 
  dplyr::slice(match(hetero_af_pos_rank$barcode, barcode)) %>% 
  tibble::column_to_rownames(var = "barcode") %>% 
  as.matrix() %>% 
  t() ->
  depth_mtx



ComplexHeatmap::Heatmap(
  matrix = depth_mtx,
  col = circlize::colorRamp2(
    breaks = range(hetero_af_pos_d$depth) %>% floor(),
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
    file = "data/PBMC_10k_v3_10x/result/04-allele-freq/mgatk_cluster_depth_heatmap-bulk.pdf",
    width = 7, 
    height = 7
  )
  ComplexHeatmap::draw(object = chm_depth)
  
  dev.off()
}

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "data/PBMC_10k_v3_10x/rda/05-mgatk-cluster.rda")
