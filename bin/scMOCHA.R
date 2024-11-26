#!/usr/bin/env Rscript
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
# library(ComplexHeatmap)
suppressPackageStartupMessages(library(ComplexHeatmap))
library(httr)
library(GetoptLong)
library(logger)
ht_opt$message <- FALSE

# src ---------------------------------------------------------------------
pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |> # nolint
  dplyr::arrange(cancer_types)


# args --------------------------------------------------------------------


# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.
#
# cell_meta_data_file <- "/home/liuc9/github/scMOCHA/06-bigdata/GSE226602/cromwell-executions/scMOCHABatch/192a6bdb-b835-4f39-a21d-9423f9c8165d/call-scMOCHA/shard-13/sub.scMOCHA/c3913f7f-efd1-4d72-9615-2463d684f359/call-gather_outputfiles/execution/GSM7080019/cell_meta_data.tsv"
# barcode_cluster_file <-"/home/liuc9/github/scMOCHA/06-bigdata/GSE226602/cromwell-executions/scMOCHABatch/192a6bdb-b835-4f39-a21d-9423f9c8165d/call-scMOCHA/shard-13/sub.scMOCHA/c3913f7f-efd1-4d72-9615-2463d684f359/call-gather_outputfiles/execution/GSM7080019/barcode_cluster.tsv"
# cell_hetero_file <-"/home/liuc9/github/scMOCHA/06-bigdata/GSE226602/cromwell-executions/scMOCHABatch/192a6bdb-b835-4f39-a21d-9423f9c8165d/call-scMOCHA/shard-13/sub.scMOCHA/c3913f7f-efd1-4d72-9615-2463d684f359/call-gather_outputfiles/execution/GSM7080019/cell.cell_heteroplasmic_df.tsv.gz"
# cell_coverage_file <-"/home/liuc9/github/scMOCHA/06-bigdata/GSE226602/cromwell-executions/scMOCHABatch/192a6bdb-b835-4f39-a21d-9423f9c8165d/call-scMOCHA/shard-13/sub.scMOCHA/c3913f7f-efd1-4d72-9615-2463d684f359/call-gather_outputfiles/execution/GSM7080019/cell.coverage.txt.gz"
# cluster_hetero_file <-"/home/liuc9/github/scMOCHA/06-bigdata/GSE226602/cromwell-executions/scMOCHABatch/192a6bdb-b835-4f39-a21d-9423f9c8165d/call-scMOCHA/shard-13/sub.scMOCHA/c3913f7f-efd1-4d72-9615-2463d684f359/call-gather_outputfiles/execution/GSM7080019/cluster.cell_heteroplasmic_df.tsv.gz"
# cluster_coverage_file <-"/home/liuc9/github/scMOCHA/06-bigdata/GSE226602/cromwell-executions/scMOCHABatch/192a6bdb-b835-4f39-a21d-9423f9c8165d/call-scMOCHA/shard-13/sub.scMOCHA/c3913f7f-efd1-4d72-9615-2463d684f359/call-gather_outputfiles/execution/GSM7080019/cluster.coverage.txt.gz"
# cell_hetero_raw_file <-"/home/liuc9/github/scMOCHA/06-bigdata/GSE226602/cromwell-executions/scMOCHABatch/192a6bdb-b835-4f39-a21d-9423f9c8165d/call-scMOCHA/shard-13/sub.scMOCHA/c3913f7f-efd1-4d72-9615-2463d684f359/call-gather_outputfiles/execution/GSM7080019/cell.cell_heteroplasmic_df_raw.tsv.gz"
# perlscript <- "/home/liuc9/github/scMOCHA/bin/get_variants_info.pl"
# jar_path <- "/scr1/users/liuc9/tools/haplogrep3"
# sqlite_path <- "/mnt/isilon/xing_lab/liuc9/refdata/mitomaster/mitomap_sqlite_20230525.sqlite3"



conda_root <- "/home/liuc9/tools/anaconda3"
conda_env <- "scmocha"
verbose <- FALSE

spec <- "
Usage: Rscript scMOCHA.R [options]

Options:
<cell_meta_data_file|meta=s> cell_meta_data.tsv
<barcode_cluster_file=s> barcode_cluster.tsv
<cell_hetero_file|ceh=s> cell.cell_heteroplasmic_df.tsv.gz
<cell_coverage_file|cec=s> cell.coverage.txt.gz
<cluster_hetero_file|clh=s> cluster.cell_heteroplasmic_df.tsv.gz
<cluster_coverage_file|clc=s> cluster.coverage.txt.gz
<cell_hetero_raw_file|chr=s> cell.cell_heteroplasmic_df_raw.tsv.gz
<perlscript=s> /home/liuc9/github/scMOCHA/bin/get_variants_info.pl
<jar_path=s> /scr1/users/liuc9/tools/haplogrep3
<sqlite_path=s> /mnt/isilon/xing_lab/liuc9/refdata/mitomaster/mitomap_sqlite_20230525.sqlite3
<conda_root=s> /home/liuc9/tools/anaconda3
<conda_env=s> scmocha
<verbose!> Print messages
"

GetoptLong.options(help_style = "two-column")
GetoptLong(spec, template_control = list(opt_width = 50))



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
    ) |>
    dplyr::filter(af > 0.05) # filter variants which AF < 0.05
}

fn_load_coverage <- function(.filename) {
  data.table::fread(
    input = .filename,
    sep = ",",
    col.names = c("pos", "barcode", "depth")
  )
}

fn_load_cluster <- function(.filename) {
  data.table::fread(
    input = .filename,
    sep = "\t",
    col.names = c("barcode", "tag", "celltype")
  ) |>
    dplyr::arrange(celltype) |>
    dplyr::mutate(celltype = factor(celltype)) |>
    dplyr::select(-tag)
}

fn_load_meta <- function(.filename) {
  data.table::fread(
    input = .filename,
    sep = "\t"
  ) |>
    dplyr::rename(
      barcode = cellbarcode
    ) |>
    dplyr::select(-orig.ident)
}

fn_af <- function(.cluster, .hetero) {
  .cluster |>
    dplyr::rename(cluster = celltype) |>
    dplyr::inner_join(
      .hetero |> tidyr::pivot_wider(
        names_from  = variant,
        values_from = af
      ),
      by = "barcode"
    )
}

fn_forplot <- function(.af, .coverage, .meta) {
  # print(.meta)
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

  .af |>
    dplyr::select(barcode, dplyr::contains(">")) |>
    tidyr::pivot_longer(
      cols = -barcode,
      names_to = "variant",
      values_to = "af"
    ) |>
    dplyr::mutate(
      pos = gsub(pattern = "([[:digit:]]*).*", "\\1", variant) |>
        as.numeric()
    ) |>
    dplyr::left_join(
      .coverage,
      by = c("barcode", "pos")
    ) |>
    tidyr::replace_na(
      replace = list(
        af = 0
      )
    ) |>
    dplyr::mutate(af = ifelse(is.na(depth), NA, af)) |>
    # dplyr::mutate(af = ifelse(depth < log2(10), -0.1, af)) |>
    dplyr::mutate(af = ifelse(depth < 10, -0.1, af)) |>
    dplyr::arrange(pos) ->
  .forplot

  .coverage |>
    dplyr::group_by(barcode) |>
    dplyr::summarise(sum_depth = sum(depth, na.rm = TRUE)) ->
  .coverage_cell

  list(
    rank = .rank,
    forplot = .forplot,
    meta = .meta,
    coverage_cell = .coverage_cell
  )
}

fn_heatmap <- function(.forplot, .cell_variants = NULL, .variant_annotation = NULL) {
  pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
    dplyr::arrange(cancer_types)
  # library(ComplexHeatmap)
  suppressPackageStartupMessages(library(ComplexHeatmap))

  .forplot$forplot |>
    dplyr::select(barcode, variant, af) |>
    tidyr::pivot_wider(
      names_from = "variant",
      values_from = af
    ) |>
    dplyr::slice(
      match(.forplot$rank$barcode, barcode)
    ) |>
    tibble::column_to_rownames(var = "barcode") |>
    as.matrix() |>
    t() ->
  .af_mtx



  tibble::tibble(
    variants = rownames(.af_mtx)
  ) ->
  .for_gcol

  .gcol <- if (is.null(.cell_variants)) {
    .for_gcol |>
      dplyr::mutate(
        cell_variants = "black"
      )
  } else {
    .for_gcol |>
      dplyr::mutate(
        cell_variants = ifelse(
          variants %in% .cell_variants,
          "black",
          "red"
        )
      )
  }


  .forplot$forplot |>
    dplyr::select(barcode, variant, depth) |>
    dplyr::mutate(
      depth = log2(depth + 1)
    ) |>
    # dplyr::arrange(pos) |>
    tidyr::pivot_wider(
      names_from = variant,
      values_from = depth
    ) |>
    dplyr::slice(match(.forplot$rank$barcode, barcode)) |>
    tibble::column_to_rownames(var = "barcode") |>
    as.matrix() |>
    t() ->
  .depth_mtx



  .forplot$rank |>
    dplyr::select(barcode, cluster) |>
    dplyr::slice(
      match(colnames(.af_mtx), barcode)
    ) |>
    dplyr::left_join(
      .forplot$meta |>
        dplyr::select(barcode, `MT%` = percent.mt),
      by = "barcode"
    ) |>
    dplyr::left_join(
      .forplot$coverage_cell |>
        dplyr::mutate(sum_depth = log10(sum_depth + 1)) |>
        dplyr::rename(`log10(Total reads)` = sum_depth),
      by = "barcode"
    ) |>
    tibble::column_to_rownames(var = "barcode") |>
    dplyr::rename(Cluster = cluster) ->
  .af_cluster


  col_clusters <- levels(.af_cluster$Cluster)
  col_colors <- pcc$color[1:length(levels(.af_cluster$Cluster))]

  names(col_colors) <- col_clusters

  chm_top <- ComplexHeatmap::HeatmapAnnotation(
    df = .af_cluster,
    # gap = unit(c(2, 2), "mm"),
    col = list(
      Cluster = col_colors,
      `MT%` = circlize::colorRamp2(
        breaks = c(2, 10),
        # colors = c("gold", "red", "black"),
        colors = c("white", "green"),
        # colors =  c("#440154FF", "#FDE725FF"),
        space = "RGB"
      ),
      `log10(Total reads)` = circlize::colorRamp2(
        # breaks = quantile(.af_cluster$`log10(Total reads)`, c(0.15, 0.75, 0.9), na.rm = T),
        breaks = quantile(.af_cluster$`log10(Total reads)`, c(0.15, 0.9), na.rm = T),
        colors = c("white", "blue"),
        space = "RGB"
      )
    ),
    which = "column"
  )

  ch_af <- if (!is.null(.variant_annotation)) {
    .df_left <- .variant_annotation |>
      dplyr::select(`Mitomap freq`, `Gnomad freq`, Haplogroup)

    .df_left |>
      dplyr::mutate(Haplogroup_col = ifelse(is.na(Haplogroup), "grey", "#3B4992FF")) |>
      dplyr::select(dplyr::contains("Haplogroup")) |>
      dplyr::filter(!is.na(Haplogroup)) |>
      dplyr::distinct() ->
    .Haplogroup

    .Haplogroup_col <- .Haplogroup$Haplogroup_col
    names(.Haplogroup_col) <- .Haplogroup$Haplogroup

    hma_left <- ComplexHeatmap::rowAnnotation(
      df = .df_left,
      col = list(
        Haplogroup = .Haplogroup_col,
        `Mitomap freq` = circlize::colorRamp2(
          breaks = c(0, 1),
          colors = c("white", "#F39B7FFF"),
          space = "RGB"
        ),
        `Gnomad freq` = circlize::colorRamp2(
          breaks = c(0, 1),
          colors = c("white", "#008280FF"),
          space = "RGB"
        )
      )
    )

    .df_right <- .variant_annotation |>
      dplyr::select(Conservation, Ntchange, Locus, Disease)

    .Ntchange <- unique(.df_right$Ntchange)
    .Ntchange_col <- rev(viridis::viridis_pal()(length(.Ntchange)))
    names(.Ntchange_col) <- .Ntchange

    hma_right <- ComplexHeatmap::rowAnnotation(
      df = .df_right,
      col = list(
        Ntchange = .Ntchange_col,
        Conservation = circlize::colorRamp2(
          breaks = c(0, 100),
          colors = c("white", "#7E6148FF"),
          space = "RGB"
        )
      )
    )

    ComplexHeatmap::Heatmap(
      matrix = .af_mtx,
      col = circlize::colorRamp2(
        breaks = c(-0.1, 0, 1),
        colors = c("lightgrey", "gold", "blue"),
        space = "RGB"
      ),
      name = "Allele Freq",
      na_col = "white",
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
      row_names_gp = gpar(
        # fontsize = 20,
        col = .gcol$cell_variants
      ),
      # column
      cluster_columns = FALSE,
      cluster_column_slices = T,
      # clustering_distance_columns = "pearson",
      # clustering_method_columns = "ward.D",
      show_column_names = FALSE,
      row_names_side = "left",
      top_annotation = chm_top,
      left_annotation = hma_left,
      right_annotation = hma_right,
      heatmap_legend_param = list(
        title = "Allele Freq",
        at = c(0, 0.5, 1),
        labels = c("0", "0.5", "1"),
        legend_direction = "vertical",
        title_gp = gpar(fontsize = 10)
      )
    )
  } else {
    ComplexHeatmap::Heatmap(
      matrix = .af_mtx,
      col = circlize::colorRamp2(
        breaks = c(-0.1, 0, 1),
        colors = c("lightgrey", "gold", "blue"),
        space = "RGB"
      ),
      name = "Allele Freq",
      na_col = "white",
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
      row_names_gp = gpar(
        # fontsize = 20,
        col = .gcol$cell_variants
      ),
      # column
      cluster_columns = FALSE,
      cluster_column_slices = T,
      # clustering_distance_columns = "pearson",
      # clustering_method_columns = "ward.D",
      show_column_names = FALSE,
      row_names_side = "left",
      top_annotation = chm_top,
      heatmap_legend_param = list(
        title = "Allele Freq",
        at = c(0, 0.5, 1),
        labels = c("0", "0.5", "1"),
        legend_direction = "vertical",
        title_gp = gpar(fontsize = 10)
      )
    )
  }


  ComplexHeatmap::Heatmap(
    matrix = .depth_mtx,
    col = circlize::colorRamp2(
      breaks = c(0, quantile(.depth_mtx, na.rm = T, probs = 0.75)),
      colors = c("white", "red"),
      # colors =  c("#440154FF", "#FDE725FF"),
      space = "RGB"
    ),
    name = "log2(Depth+1)",
    na_col = "white",
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
    row_names_gp = gpar(
      # fontsize = 20,
      col = .gcol$cell_variants
    ),
    # column
    cluster_columns = FALSE,
    cluster_column_slices = T,
    # clustering_distance_columns = "pearson",
    # clustering_method_columns = "ward.D",
    show_column_names = FALSE,
    row_names_side = "left",
    top_annotation = chm_top
  ) ->
  ch_depth

  list(
    ch_af = ch_af,
    ch_depth = ch_depth
  )
}

fn_plot_cell_violin <- function(.forplot, .cell_anno) {
  pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
    dplyr::arrange(cancer_types)

  .forplot$forplot |>
    dplyr::filter(af > 0) |>
    dplyr::filter(variant %in% .cell_anno$variant) |>
    dplyr::left_join(
      .forplot$rank |> dplyr::select(-s_af),
      by = "barcode"
    ) ->
  .theforplot

  .theforplot |>
    dplyr::group_by(variant) |>
    dplyr::summarise(maf = mean(af, na.rm = T)) |>
    dplyr::arrange(-maf) ->
  .sort_variant

  .theforplot |>
    dplyr::group_by(cluster, variant) |>
    dplyr::summarise(mean_cluster_variant_af = mean(af, na.rm = T)) |>
    dplyr::ungroup() ->
  .cluster_variant_af

  .cell_anno |>
    dplyr::filter(variant %in% .sort_variant$variant) |>
    dplyr::mutate(fill = ifelse(!is.na(Haplogroup), "#3B0049", "white")) |>
    dplyr::mutate(color = ifelse(!is.na(Haplogroup), "white", "black")) |>
    dplyr::mutate(
      variant = factor(variant, .sort_variant$variant)
    ) |>
    dplyr::arrange(variant) ->
  .haplo_variant

  .theforplot |>
    dplyr::left_join(
      .cluster_variant_af,
      by = c("cluster", "variant")
    ) |>
    dplyr::mutate(
      variant = factor(variant, .sort_variant$variant |> unique())
    ) |>
    dplyr::arrange(variant) ->
  .haplo_forplot


  library(ggh4x)
  library(ggbeeswarm)
  .haplo_forplot |>
    ggplot(aes(x = cluster, y = af)) +
    geom_violin(
      aes(fill = mean_cluster_variant_af),
      alpha = 0.5,
      size = 1,
      color = NA
    ) +
    ggbeeswarm::geom_quasirandom(
      # shape = 21,
      aes(color = af),
      size = 1,
      dodge.width = .75,
      alpha = .5,
      show.legend = T,
    ) +
    ggh4x::facet_wrap2(
      ~variant,
      ncol = 12,
      strip = ggh4x::strip_themed(
        background_x = elem_list_rect(
          fill = .haplo_variant$fill
        ),
        text_x = elem_list_text(
          colour = .haplo_variant$color,
          face = c("bold")
        ),
        by_layer_x = FALSE,
      )
    ) +
    scale_fill_gradient2(
      name = "AF",
      low = "white",
      mid = "red",
      high = "#3B0049",
      midpoint = 0.5
    ) +
    scale_color_gradient2(
      name = "AF",
      low = "white",
      mid = "red",
      high = "#3B0049",
      midpoint = 0.5
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(
        color = "black",
      ),
      legend.position = "none ",
      plot.title = element_text(
        size = 16,
        hjust = 0.5
      ),
      axis.line = element_line(
        color = "black"
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        # color = pcc$color
      )
    ) ->
  p
  p
}




# cell cluster ------------------------------------------------------------

log_info("load cluster_umap")
cluster_umap <- fn_load_cluster(
  .filename = barcode_cluster_file
)

log_info("load metadata")
metadata <- fn_load_meta(
  .filename = cell_meta_data_file
)
# Cell allele -------------------------------------------------------------

cell_hetero <- fn_load_hetero(
  .filename = cell_hetero_file
)

cell_coverage <- fn_load_coverage(
  .filename = cell_coverage_file
)

cell_cluster_af <- fn_af(
  .cluster = cluster_umap,
  .hetero = cell_hetero
)

cell_cluster_forplot <- fn_forplot(
  .af = cell_cluster_af,
  .coverage = cell_coverage,
  .meta = metadata
)

log_info("fn_heatmap")
# print(cell_cluster_forplot)

ch_af_depth <- fn_heatmap(
  .forplot = cell_cluster_forplot,
  .cell_variants = NULL,
  .variant_annotation = NULL
)



{
  pdf(
    file = "cell_af_heatmap.pdf",
    width = 14,
    height = 15
  )
  ComplexHeatmap::draw(object = ch_af_depth$ch_af)
  dev.off()

  pdf(
    file = "cell_depth_heatmap.pdf",
    width = 14,
    height = 15
  )
  ComplexHeatmap::draw(object = ch_af_depth$ch_depth)
  dev.off()
  log_success("save image")
}


# cluster allele-----------------------------------------------------------------


cluster_hetero <- fn_load_hetero(
  .filename = cluster_hetero_file
) |>
  dplyr::mutate(cluster = barcode) |>
  dplyr::mutate(cluster = factor(cluster)) |>
  dplyr::left_join(
    cluster_umap |>
      dplyr::mutate(cluster = celltype) |>
      dplyr::mutate(cluster = factor(cluster)) |>
      dplyr::select(cluster, celltype) |>
      dplyr::distinct(),
    by = "cluster"
  ) |>
  dplyr::select(-cluster) |>
  dplyr::rename(cluster = celltype)

cluster_coverage <- fn_load_coverage(
  .filename = cluster_coverage_file
)


cluster_cluster_af <-
  cluster_hetero |> tidyr::pivot_wider(
    names_from  = variant,
    values_from = af
  )

cluster_cluster_forplot <- fn_forplot(
  .af = cluster_cluster_af,
  .coverage = cluster_coverage,
  .meta = metadata
)


cluster_ch_af_depth <- fn_heatmap(
  .forplot = cluster_cluster_forplot
)


{
  pdf(
    file = "cluster_af_heatmap.pdf",
    width = 7,
    height = 15
  )
  ComplexHeatmap::draw(object = cluster_ch_af_depth$ch_af)
  dev.off()

  pdf(
    file = "cluster_depth_heatmap.pdf",
    width = 7,
    height = 15
  )
  ComplexHeatmap::draw(object = cluster_ch_af_depth$ch_depth)
  dev.off()
}

venn_cell_cluster <- ggvenn::ggvenn(
  data = list(
    Cell = unique(cell_hetero$variant),
    Cluster = unique(cluster_hetero$variant)
  ),
  fill_color = ggsci::pal_aaas()(2)
)

ggsave(
  filename = "venn_cell_cluster.pdf",
  plot = venn_cell_cluster,
  device = "pdf",
  width = 7,
  height = 5
)


# Cluster cell allele -----------------------------------------------------

cell_hetero_raw <- fn_load_hetero(
  .filename = cell_hetero_raw_file
) |>
  dplyr::filter(
    variant %in% cluster_hetero$variant # only keep the variants in cluster_hetero
    # variant %in% c(cell_hetero$variant, cluster_hetero$variant)
  )

cell_raw_cluster_af <- cluster_umap |>
  dplyr::left_join(cell_hetero_raw, by = "barcode") |>
  dplyr::rename(cluster = celltype) |>
  tidyr::pivot_wider(
    names_from = variant,
    values_from = af
  )

cell_raw_cluster_forplot <- fn_forplot(
  .af = cell_raw_cluster_af,
  .coverage = cell_coverage,
  .meta = metadata
)


# Variant annotation ------------------------------------------------------


cell_raw_cluster_forplot$forplot |>
  dplyr::filter(!is.na(depth)) |>
  # dplyr::select(barcode, pos, variant) |>
  dplyr::select(pos, variant) |>
  dplyr::distinct() |>
  dplyr::mutate(variant = gsub(
    pattern = "[0-9]*",
    replacement = "",
    x = variant
  )) |>
  tidyr::separate(
    col = variant,
    into = c("ref", "var")
  ) |>
  # dplyr::rename(sample = barcode) |>
  dplyr::mutate(sample = "sample1") |>
  dplyr::select(
    sample = sample,
    pos = pos,
    ref = ref,
    var = var
  ) |>
  dplyr::mutate(
    v = glue::glue("{pos}{ref}>{var}")
  ) |>
  dplyr::select(sample, v) |>
  tibble::rowid_to_column() |>
  tidyr::pivot_wider(
    names_from = rowid,
    values_from = v
  ) ->
cell_variants

readr::write_delim(
  x = cell_variants,
  file = "cell_snvlist.tsv",
  delim = " ",
  col_names = FALSE
)


cmd <- "source {conda_root}/etc/profile.d/conda.sh; conda activate {conda_env}; perl {perlscript} {file.path(jar_path, 'haplogrep3.jar')} {sqlite_path} cell_snvlist.tsv > cell_variant_annotation.tsv" |> glue::glue()
log_debug(cmd)
system(command = cmd)

if (file.exists("cell_variant_annotation.tsv")) {
  cell_anno <- readr::read_tsv("cell_variant_annotation.tsv") |>
    dplyr::mutate(variant = glue::glue("{Position}{Ref}>{Alt}"))

  writexl::write_xlsx(
    x = cell_anno,
    path = "cell_variant_annotation.xlsx"
  )


  variant_annotation <- cell_anno |>
    dplyr::mutate(
      variant = glue::glue("{Position}{Ref}>{Alt}")
    ) |>
    dplyr::mutate(
      Status = ifelse(
        !is.na(Status),
        "Reported",
        Status
      )
    ) |>
    dplyr::select(
      variant, ntchange,
      calc_locus = Locus,
      Haplogroup,
      Verbose_haplogroup,
      Disease,
      Status,
      Conservation,
      mito_freq = `Mitomap Frequency`,
      gnomad_freq = `Gnomad Frequency`
    ) |>
    dplyr::mutate(
      calc_locus = gsub(
        pattern = "<br>.*",
        replace = "",
        x = calc_locus
      )
    ) |>
    dplyr::mutate(
      Conservation = gsub(
        pattern = "%",
        replacement = "",
        x = Conservation
      )
    ) |>
    dplyr::mutate(
      Disease = stringr::str_wrap(
        stringr::str_to_sentence(string = Disease),
        width = 30
      )
    ) |>
    dplyr::mutate(Conservation = as.numeric(Conservation)) |>
    dplyr::mutate(
      mito_freq = mito_freq / 100,
      gnomad_freq = gnomad_freq / 100
    ) |>
    dplyr::select(
      Ntchange = ntchange,
      Locus = calc_locus,
      Haplogroup = Verbose_haplogroup,
      Disease = Disease,
      Status,
      Conservation,
      `Mitomap freq` = mito_freq,
      `Gnomad freq` = gnomad_freq
    )
} else {
  variant_annotation <- NULL
}


cell_raw_ch_af_depth <- fn_heatmap(
  .forplot = cell_raw_cluster_forplot,
  # .cell_variants = cell_cluster_forplot$forplot$variant,
  .cell_variants = unique(cell_hetero$variant),
  .variant_annotation = variant_annotation
)

{
  pdf(
    file = "cluster_cell_af_heatmap.pdf",
    width = 25,
    height = 15
  )
  ComplexHeatmap::draw(object = cell_raw_ch_af_depth$ch_af)
  dev.off()

  pdf(
    file = "cluster_cell_depth_heatmap.pdf",
    width = 14,
    height = 15
  )
  ComplexHeatmap::draw(object = cell_raw_ch_af_depth$ch_depth)
  dev.off()
}

# violin plot -------------------------------------------------------------



fn_plot_cell_violin(
  .forplot = cell_raw_cluster_forplot,
  .cell_anno = cell_anno
) -> p_violin

{
  pdf(
    file = "cluster_cell_violin.pdf",
    width = 24,
    height = 12
  )
  print(p_violin)
  dev.off()
}



# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------

save.image(file = "scMOCHA.rda")
# load(file = "/scr1/users/liuc9/tmp/mito/flu2-a/cromwell-executions/SCMTAH/0138fcd0-c384-42c2-8704-6647767610d2/call-plot_scmtah/execution/scmtah.rda")
