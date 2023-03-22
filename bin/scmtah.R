# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed Mar 22 15:04:30 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

# src ---------------------------------------------------------------------

args <- commandArgs(TRUE)

sc_dir <- "/home/liuc9/scratch/mitochondrial/testdata/1_old_donor_pbmc/flu2/Flu2/outs"
cluster_dir <- "/home/liuc9/scratch/mitochondrial/testdata/1_old_donor_pbmc/flu2/Flu2/outs" 

# hetero_file <- args[1]
# coverage_file <- args[2]
# cluster_umap_file <- args[3]


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------

hetero <- data.table::fread(input = file.path(
  sc_dir,
  "mgatk_out/final",
  "sc.cell_heteroplasmic_df.tsv.gz"
)) |> 
  dplyr::rename(barcode = "V1") |> 
  tidyr::pivot_longer(
    cols = -barcode,
    names_to = "variant",
    values_to = "af"
  )

coverage <- data.table::fread(
  input = file.path(
    sc_dir,
    "mgatk_out/final",
    "sc.coverage.txt.gz"
  ),
  sep = ",",
  col.names = c("pos", "barcode", "depth")
) |> 
  dplyr::mutate(depth = log2(depth + 1))

coverage_wider <- coverage |> 
  tidyr::pivot_wider(
    names_from = pos,
    values_from = depth
  )

cluster_umap <- data.table::fread(
  input = file.path(
    sc_dir,
    "cluster_umap.tsv"
  ),
  sep = "\t"
) |> 
  dplyr::mutate(cluster = factor(cluster))

# body --------------------------------------------------------------------

cluster_umap %>%
  dplyr::left_join(
    hetero |> tidyr::pivot_wider(names_from  = variant, values_from = af), 
    by = "barcode"
  ) ->
  cell_cluster_af


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


coverage |> 
  dplyr::filter(barcode %in% cell_cluster_af_pos$barcode) %>%
  dplyr::filter(pos %in% cell_cluster_af_pos$pos) ->
  coverage_pos


cell_cluster_af_pos %>%
  dplyr::left_join(
    coverage_pos,
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
  dplyr::select(barcode, cluster) %>%
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



cell_cluster_af_pos %>%
  dplyr::left_join(
    coverage_pos,
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
    breaks = c(0, 4),
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




# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------