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
library(ComplexHeatmap)
pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv")

# src ---------------------------------------------------------------------

args <- commandArgs(TRUE)

sc_dir <- "/home/liuc9/scratch/mitochondrial/testdata/1_old_donor_pbmc/flu2/Flu2/outs"
# cluster_dir <- "/home/liuc9/scratch/mitochondrial/testdata/1_old_donor_pbmc/flu2/Flu2/outs" 

# hetero_file <- args[1]
# coverage_file <- args[2]
# cluster_umap_file <- args[3]


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------

fn_load_hetero <- function(.filename) {
  # .filename <- file.path(
  #   sc_dir,
  #   "mgatk_out/final",
  #   "sc.cell_heteroplasmic_df.tsv.gz"
  # )
  
  data.table::fread(input = .filename) |> 
    dplyr::rename(barcode = "V1") |> 
    tidyr::pivot_longer(
      cols = -barcode,
      names_to = "variant",
      values_to = "af"
    )
}

fn_load_coverage <- function(.filename) {
  
  data.table::fread(
    input = .filename,
    sep = ",",
    col.names = c("pos", "barcode", "depth")
  ) |> 
    dplyr::mutate(depth = log2(depth + 1))
}

fn_load_cluster <- function(.filename) {
  data.table::fread(
    input = .filename,
    sep = "\t"
  ) |> 
    dplyr::mutate(cluster = factor(cluster))
}

fn_af <- function(.cluster, .hetero) {
  .cluster %>%
    dplyr::left_join(
      .hetero |> tidyr::pivot_wider(
        names_from  = variant, 
        values_from = af
      ), 
      by = "barcode"
    ) 
}

fn_forplot <- function(.af, .coverage) {
  .af |> 
    dplyr::select(barcode, cluster, dplyr::contains(">")) |> 
    tidyr::pivot_longer(
      cols = -c(barcode, cluster),
      names_to = "variant",
      values_to = "af"
    ) |> 
    dplyr::group_by(barcode, cluster) |> 
    dplyr::summarise(s_af = sum(af, na.rm = T)) |> 
    dplyr::ungroup() |> 
    dplyr::arrange(cluster, -s_af) ->
    .rank
  
  .af %>%
    dplyr::select(barcode, dplyr::contains(">")) %>%
    tidyr::pivot_longer(
      cols = -barcode,
      names_to = "variant",
      values_to = "af"
    ) %>%
    dplyr::mutate(
      pos = gsub(pattern = "([[:digit:]]*).*", "\\1", variant) %>%
        as.numeric()
    ) |> 
    dplyr::left_join(
      .coverage,
      by = c("barcode", "pos")
    ) %>%
    tidyr::replace_na(
      replace = list(
        af = 0
      )
    ) %>%
    dplyr::mutate(af = ifelse(is.na(depth), NA, af)) %>%
    dplyr::arrange(pos) ->
    .forplot
  
  list(
    rank = .rank,
    forplot = .forplot
  )
}

fn_heatmap <- function(.forplot) {
  
  .forplot$forplot |> 
    dplyr::select(barcode, variant, af) %>%
    tidyr::pivot_wider(
      names_from = "variant",
      values_from = af
    ) %>%
    dplyr::slice(
      match(.forplot$rank$barcode, barcode)
    ) %>%
    tibble::column_to_rownames(var = "barcode") %>%
    as.matrix() %>%
    t() ->
    .af_mtx
  
  .forplot$forplot |> 
    dplyr::select(barcode, pos, depth) %>%
    dplyr::arrange(pos) %>%
    tidyr::pivot_wider(
      names_from = pos,
      values_from = depth
    ) %>%
    dplyr::slice(match(.forplot$rank$barcode, barcode)) %>%
    tibble::column_to_rownames(var = "barcode") %>%
    as.matrix() %>%
    t() ->
    .depth_mtx
  
  
  .forplot$rank %>%
    dplyr::select(barcode, cluster) %>%
    dplyr::slice(
      match(colnames(.af_mtx), barcode)
    ) %>%
    tibble::column_to_rownames(var = "barcode") %>%
    dplyr::select(Cluster = cluster) ->
    .af_cluster
  
  
  col_clusters <- levels(.af_cluster$Cluster) %>% as.numeric()
  col_colors <- pcc$color[1:length(levels(.af_cluster$Cluster))]
  
  names(col_colors) <- col_clusters
  
  chm_top <- ComplexHeatmap::HeatmapAnnotation(
    df = .af_cluster,
    gap = unit(c(2,2), "mm"),
    col = list(Cluster = col_colors),
    which = "column"
  )
  
  
  ComplexHeatmap::Heatmap(
    matrix = .af_mtx,
    col = circlize::colorRamp2(
      breaks = c(0, 0.98, 1),
      colors = c("white", "#440154FF", "#FDE725FF"),
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
    ch_af
  
  
  ComplexHeatmap::Heatmap(
    matrix = .depth_mtx,
    col = circlize::colorRamp2(
      breaks = c(0, quantile(.depth_mtx, na.rm = T, probs = 0.75)),
      colors = c("white", "red"),
      # colors =  c("#440154FF", "#FDE725FF"),
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
    ch_depth
  
  list(
    ch_af = ch_af,
    ch_depth = ch_depth
  )
}

# cell cluster ------------------------------------------------------------

cluster_umap <- fn_load_cluster(
  .filename = file.path(
    sc_dir,
    "cluster_umap.tsv"
  )
)

# Cell allele -------------------------------------------------------------
hetero <- fn_load_hetero(
  .filename = file.path(
    sc_dir,
    "mgatk_out/final",
    "sc.cell_heteroplasmic_df.tsv.gz"
  )
)

coverage <- fn_load_coverage(
  .filename = file.path(
    sc_dir,
    "mgatk_out/final",
    "sc.coverage.txt.gz"
  ) 
)



cell_cluster_af <- fn_af(.cluster = cluster_umap, .hetero = hetero)

cell_cluster_forplot <- fn_forplot(.af = cell_cluster_af, .coverage = coverage)

ch_af_depth <- fn_heatmap(.forplot = cell_cluster_forplot)


ch_af_depth$ch_af
ch_af_depth$ch_depth

{
  pdf(
    file = "mgatk_cell_al_heatmap.pdf",
    width = 14, 
    height = 7
  )
  ComplexHeatmap::draw(object = ch_af_depth$ch_af)
  dev.off()

  pdf(
    file = "mgatk_cell_depth_heatmap.pdf",
    width = 14, 
    height = 7
  )
  ComplexHeatmap::draw(object = ch_af_depth$ch_depth)
  dev.off()
}


# cluster allele-----------------------------------------------------------------


hetero_cluster <- fn_load_hetero(
  .filename = file.path(
    sc_dir,
    "mgatk_cluster/final",
    "mgatk_cluster.cell_heteroplasmic_df.tsv.gz"
  )
) |> 
  dplyr::mutate(
    barcode = gsub(
      pattern = "cluster|-1",
      replacement = "",
      x = barcode
    )
  ) %>% 
  dplyr::mutate(barcode = as.integer(barcode) -1) |> 
  dplyr::mutate(cluster = barcode) |> 
  dplyr::mutate(cluster = factor(cluster))

coverage_cluster <- fn_load_coverage(
  .filename = file.path(
    sc_dir,
    "mgatk_cluster/final",
    "mgatk_cluster.coverage.txt.gz"
  ) 
) |> 
  dplyr::mutate(
    barcode = gsub(
      pattern = "cluster|-1",
      replacement = "",
      x = barcode
    )
  ) %>% 
  dplyr::mutate(barcode = as.integer(barcode) -1)


cluster_cluster_af <- 
  hetero_cluster |> tidyr::pivot_wider(
      names_from  = variant, 
      values_from = af
  )

cluster_cluster_forplot <- fn_forplot(
  .af = cluster_cluster_af, 
  .coverage = coverage_cluster
  )


cluster_ch_af_depth <- fn_heatmap(.forplot = cluster_cluster_forplot)

{
  pdf(
    file = "mgatk_cluster_al_heatmap.pdf",
    width = 7, 
    height = 7
  )
  ComplexHeatmap::draw(object = cluster_ch_af_depth$ch_af)
  dev.off()
  
  pdf(
    file = "mgatk_cluster_depth_heatmap.pdf",
    width = 7, 
    height = 7
  )
  ComplexHeatmap::draw(object = cluster_ch_af_depth$ch_depth)
  dev.off()
}


# Cluster allele ----------------------------------------------------------



# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------