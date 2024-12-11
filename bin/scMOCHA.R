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

fn_heatmap <- function(.forplot, .cell_variants = NULL, .variant_annotation = NULL, col_option = "turbo", show_column_title = FALSE) {
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




  # col_start = 0
  # col_end = 1
  # n_break = 100
  sort(unique(as.numeric(.af_mtx))) -> .seq_af
  # .seq_af |> hist()
  # .seq_af_sub <- c(
  #   -0.1, 0, 0.05, 0.1, 0.15, 0.2,
  #   .seq_af[(.seq_af > 0.2 & .seq_af < 0.8)],
  #   0.8, 0.85, 0.9, 0.95, 1
  # )
  # .seq_af_sub |> hist()
  if (show_column_title) {
    .column_title <- paste0("palette = '", col_option, "'")
  } else {
    .column_title <- NULL
  }

  col_pick = c(
    "grey",
    viridis::viridis_pal(
      alpha = 1,
      begin = 0,
      end = 1,
      direction = 1,
      option = col_option
    )(
      # n_break
      length(.seq_af) - 1
    )
  )
  # (a) change the “0” value from dark blue to dark green;
  # col_pick |> color()
  #  (b) reverse the color, use red as “0” and violet as “1”
  col_fun = circlize::colorRamp2(
    # seq(col_start,
    #   col_end,
    #   length.out = n_break
    # ),
    .seq_af,
    col_pick
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
    # .Ntchange_col <- rev(viridis::viridis_pal()(length(.Ntchange)))
    .Ntchange_col <- c("#FDE725FF", "#440154FF")
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
      # col = circlize::colorRamp2(
      #   breaks = c(-0.1, 0, 1),
      #   colors = c("lightgrey", "gold", "blue"),
      #   space = "RGB"
      # ),
      col = col_fun,
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
      # column_title = paste0("palette = '", col_option, "'"),
      column_title = .column_title,
      column_title_gp = gpar(fontsize = 40),
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
    # col_fun = circlize::colorRamp2(c(0, 0.5, 1), hcl_palette = "turbo")

    ComplexHeatmap::Heatmap(
      matrix = .af_mtx,
      # col = circlize::colorRamp2(
      #   breaks = c(-0.1, 0, 1),
      #   colors = c("lightgrey", "gold", "blue"),
      #   space = "RGB"
      # ),
      col = col_fun,
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
      # column_title = paste0("palette = '", col_option, "'"),
      column_title = .column_title,
      column_title_gp = gpar(fontsize = 40),
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

fn_plot_cell_violin <- function(.forplot, .cell_anno, .sel_variants = NULL) {
  pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
    dplyr::arrange(cancer_types)

  if (!is.null(.sel_variants)) {
    .cell_anno <- .cell_anno |>
      dplyr::filter(variant %in% .sel_variants)
  }

  .forplot$forplot |>
    dplyr::left_join(
      .forplot$rank |> dplyr::select(-s_af),
      by = "barcode"
    ) |>
    dplyr::filter(variant %in% .cell_anno$variant) |>
    dplyr::mutate(
      af = ifelse(af == 0, NA_real_, af)
    ) |>
    dplyr::mutate(
      af = ifelse(depth < 10, NA_real_, af)
    ) |>
    dplyr::mutate(
      depth_log2 = log2(depth + 1)
    ) ->
  .forplot_cluster_cell_variant

  .forplot_cluster_cell_variant |>
    dplyr::filter(af > 0) ->
  .theforplot

  .theforplot |>
    dplyr::select(variant, pos) |>
    dplyr::distinct() |>
    dplyr::arrange(pos) ->
  .sort_variant

  .forplot_cluster_cell_variant |>
    dplyr::group_by(cluster, variant) |>
    dplyr::summarise(
      mean_cluster_variant_af = mean(af, na.rm = T),
      sum_cluster_variant_depth = sum(depth, na.rm = T),
      max_cluster_variant_depth = max(depth, na.rm = T)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      sum_cluster_variant_depth_log2 = log2(sum_cluster_variant_depth + 1)
    ) |>
    dplyr::filter(variant %in% .sort_variant$variant) ->
  .cluster_variant_af

  .c <- unique(.cluster_variant_af$cluster)
  .v <- unique(.cluster_variant_af$variant)

  tibble::tibble(
    cluster = rep(.c, each = length(.v)),
    variant = rep(.v, length(.c))
  ) |>
    dplyr::left_join(
      .cluster_variant_af,
      by = c("cluster", "variant")
    ) |>
    dplyr::mutate(
      rect_type = ifelse(max_cluster_variant_depth >= 10, "white", "grey")
    ) |>
    dplyr::mutate(
      variant = factor(variant, .sort_variant$variant)
    ) ->
  .no_depth

  # .theforplot |>
  #   dplyr::group_by(variant) |>
  #   dplyr::summarise(maf = mean(af, na.rm = T)) |>
  #   dplyr::arrange(-maf) ->
  # .sort_variant

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
    dplyr::inner_join(
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
  library(ggnewscale)

  .haplo_forplot |>
    ggplot(aes(x = cluster)) +
    ggh4x::facet_wrap2(
      ~variant,
      ncol = 12,
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = .haplo_variant$fill
        ),
        text_x = ggh4x::elem_list_text(
          colour = .haplo_variant$color,
          face = c("bold")
        ),
        by_layer_x = FALSE,
      )
    ) +
    geom_rect(
      data = .no_depth,
      aes(
        xmin = as.numeric(cluster) - 0.5,
        xmax = as.numeric(cluster) + 0.5,
        ymin = -Inf,
        ymax = Inf,
        fill = rect_type
      ),
      alpha = 0.2,
      color = NA,
      show.legend = FALSE
    ) +
    scale_fill_identity() +
    new_scale_fill() +
    geom_violin(
      aes(
        y = af,
        fill = mean_cluster_variant_af
      ),
      alpha = 0.5,
      size = 1,
      color = NA,
      show.legend = FALSE
    ) +
    scale_fill_gradient2(
      name = "AF",
      low = "white",
      mid = "red",
      high = "#3B0049",
      midpoint = 0.5,
    ) +
    ggbeeswarm::geom_quasirandom(
      aes(
        y = af,
        color = af
      ),
      size = 1,
      dodge.width = .75,
      alpha = .5,
    ) +
    scale_color_gradient2(
      name = "AF",
      low = "white",
      mid = "red",
      high = "#3B0049",
      midpoint = 0.5,
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(
        color = "black",
      ),
      # legend.position = "none ",
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
  p_af

  .forplot_cluster_cell_variant |>
    dplyr::filter(variant %in% .sort_variant$variant) |>
    dplyr::left_join(
      .cluster_variant_af,
      by = c("cluster", "variant")
    ) |>
    dplyr::mutate(
      variant = factor(variant, .sort_variant$variant)
    ) |>
    # .haplo_forplot |>
    ggplot(aes(x = cluster)) +
    ggh4x::facet_wrap2(
      ~variant,
      ncol = 12,
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = .haplo_variant$fill
        ),
        text_x = ggh4x::elem_list_text(
          colour = .haplo_variant$color,
          face = c("bold")
        ),
        by_layer_x = FALSE,
      )
    ) +
    geom_rect(
      data = .no_depth,
      aes(
        xmin = as.numeric(cluster) - 0.5,
        xmax = as.numeric(cluster) + 0.5,
        ymin = -Inf,
        ymax = Inf,
        fill = rect_type
      ),
      alpha = 0.2,
      color = NA,
      show.legend = FALSE
    ) +
    scale_fill_identity() +
    new_scale_fill() +
    geom_violin(
      aes(
        y = depth_log2,
        fill = sum_cluster_variant_depth_log2
      ),
      alpha = 0.5,
      size = 1,
      color = NA,
      show.legend = FALSE
    ) +
    scale_fill_gradient2(
      name = "log2(Depth+1)",
      low = "white",
      mid = "red",
      high = "#3B0049",
      midpoint = 0.5,
    ) +
    ggbeeswarm::geom_quasirandom(
      # shape = 21,
      aes(
        y = depth_log2,
        color = depth_log2
      ),
      size = 1,
      dodge.width = .75,
      alpha = .5,
    ) +
    scale_color_gradient2(
      name = "log2(Depth+1)",
      low = "white",
      mid = "red",
      high = "#3B0049",
      midpoint = 0.5,
    ) +
    geom_hline(
      yintercept = log2(10),
      linetype = 21,
      color = "black"
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(
        color = "black",
      ),
      # legend.position = "none ",
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
  p_depth


  list(
    p_af = p_af,
    p_depth = p_depth,
    haplo_variant = .haplo_variant,
    haplo_forplot = .haplo_forplot
  )
}

fn_somatic_variant <- function(.haplo_variant, .haplo_violin, .n_cells = 10, .high_af = 0.95) {
  # .haplo_variant <- srr_out_cell_stats$haplo_variant[[27]]
  # .haplo_violin <- srr_out_cell_stats$haplo_violin[[27]]

  # 1. filter by haplogrep marker variant
  .haplo_variant |>
    dplyr::filter(fill != "white") |>
    dplyr::pull(variant) ->
  .v_haplo

  # 2. filter by n_cells
  .n_cells <- 10
  .haplo_violin |>
    dplyr::count(variant) |>
    dplyr::filter(n < .n_cells) |>
    dplyr::pull(variant) ->
  .v_n_cells

  # 3. tRNA p9 and RNA editing position
  .editing_pos <- c(
    585, 1610, 3238, 4271, 5520, 7526, 8303, # tRNA p9
    9999, 10413, 12146, 12274, 14734, 15896, # tRNA p9
    295, 2617, 13710 # RNA editing
  )

  .haplo_variant |>
    dplyr::filter(Position %in% .editing_pos) |>
    dplyr::pull(variant) ->
  .v_editing

  # 4. high af in 95% of cells
  .haplo_violin |>
    dplyr::group_by(variant) |>
    dplyr::summarise(
      afm = mean(af, na.rm = T)
    ) |>
    dplyr::filter(afm > .high_af) |>
    dplyr::pull(variant) ->
  .v_high_af

  # somatic variant
  .haplo_variant |>
    dplyr::filter(!variant %in% c(.v_haplo, .v_n_cells, .v_editing, .v_high_af)) |>
    dplyr::pull(variant) ->
  .v_somatic

  list(
    haplo = .v_haplo,
    n_cells = .v_n_cells,
    editing = .v_editing,
    somatic = .v_somatic,
    high_af = .v_high_af
  )
}

fn_plot_mtdna <- function() {
  mt_exons_df <- "/home/liuc9/github/scMOCHA/fasta/mt_exons.df.rds.gz"


  gtf_gene_df <-
    readr::read_rds(
      file = mt_exons_df
    )
  library(gggenes)
  ggplot(gtf_gene_df, aes(xmin = start, xmax = end, y = seqnames)) +
    # geom_gene_arrow() +
    geom_gene_arrow(
      aes(
        fill = gene_biotype
      ),
      arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")
    ) +
    scale_fill_brewer(
      palette = "Set1",
      name = "Gene type",
      labels = c("MT rRNA", "MT tRNA", "Protein coding")
    ) +
    ggrepel::geom_text_repel(
      aes(x = (start + end) / 2, label = gene_name, color = gene_biotype),
      # fill = "white",
      # nudge_x =1,
      # nudge_y = -0.1,
      size = 3,
      show.legend = F,
      max.overlaps = Inf,
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_x_continuous(
      limits = c(0, 17000),
      breaks = seq(0, 17000, 1000),
      expand = expansion(mult = c(0, 0.03)),
    ) +
    scale_y_discrete(
      expand = expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    # theme_genes() +
    theme(
      legend.position = "bottom",
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      # axis.text.x = element_text(size = 14),
      # legend.text = element_text(size = 14),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.text.x = element_text(
        vjust = -1,
      ),
    ) ->
  pg
  pg
}

fn_plot_coverage <- function(.cluster_coverage) {
  pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
    dplyr::arrange(cancer_types)

  .cluster_coverage |>
    ggplot(aes(x = pos, y = depth, fill = barcode)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_x_continuous(
      expand = expansion(mult = c(0.01, 0)),
      limits = c(1, 17000),
      breaks = seq(0, 17000, 1000),
      labels = seq(0, 17000, 1000)
    ) +
    scale_y_continuous(
      expand = c(0.01, 0),
      # limits = c(0, 520000),
      label = scales::label_number()
    ) +
    scale_fill_manual(
      name = "Cell type",
      values = pcc$color
    ) +
    ggh4x::facet_wrap2(
      ~barcode,
      ncol = 1,
      strip.position = "right",
      strip = ggh4x::strip_themed(
        background_y = ggh4x::elem_list_rect(
          fill = pcc$color
        ),
        text_y = ggh4x::elem_list_text(
          colour = "white",
          face = c("bold")
        ),
        by_layer_y = FALSE,
      )
    ) +
    theme(
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line.y.left = element_line(color = "black"),
      # axis.line.x.bottom = element_line(color = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = c(0.8, 0.5),
      legend.key = element_blank(),
      axis.title.y = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      legend.text = element_text(
        size = 14,
        color = "black"
      ),
      legend.title = element_text(
        size = 16,
        colour = "black"
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        size = 8,
        color = "black",
        face = "bold"
      )
    ) +
    labs(y = "Depth") ->
  p_mt_depth_celltype

  .cluster_coverage |>
    dplyr::group_by(pos) |>
    dplyr::summarise(depth = sum(depth, na.rm = T)) |>
    ggplot(aes(x = pos, y = depth)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_x_continuous(
      expand = expansion(mult = c(0.01, 0)),
      limits = c(1, 17000),
      breaks = seq(0, 17000, 1000),
      labels = seq(0, 17000, 1000)
    ) +
    scale_y_continuous(
      expand = c(0.01, 0),
      # limits = c(0, 520000),
      label = scales::label_number()
    ) +
    theme(
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line.y.left = element_line(color = "black"),
      # axis.line.x.bottom = element_line(color = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = c(0.8, 0.5),
      legend.key = element_blank(),
      axis.title.y = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      legend.text = element_text(
        size = 14,
        color = "black"
      ),
      legend.title = element_text(
        size = 16,
        colour = "black"
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        size = 8,
        color = "black",
        face = "bold"
      )
    ) +
    labs(y = "Depth") ->
  p_mt_depth_allcell

  list(
    p_mt_depth_celltype = p_mt_depth_celltype,
    p_mt_depth_allcell = p_mt_depth_allcell
  )
}

fn_plot_hotspots <- function(.forplot, .cell_anno, .sel_variants = NULL) {
  pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
    dplyr::arrange(cancer_types)

  if (!is.null(.sel_variants)) {
    .cell_anno <- .cell_anno |>
      dplyr::filter(variant %in% .sel_variants)
  }

  .forplot$forplot |>
    dplyr::left_join(
      .forplot$rank |> dplyr::select(-s_af),
      by = "barcode"
    ) |>
    dplyr::filter(variant %in% .cell_anno$variant) |>
    dplyr::mutate(
      af = ifelse(af == 0, NA_real_, af)
    ) |>
    dplyr::mutate(
      af = ifelse(depth < 10, NA_real_, af)
    ) |>
    dplyr::mutate(
      depth_log2 = log2(depth + 1)
    ) ->
  .forplot_cluster_cell_variant

  .forplot_cluster_cell_variant |>
    dplyr::filter(af > 0) ->
  .theforplot

  .theforplot |>
    dplyr::select(variant, pos) |>
    dplyr::distinct() |>
    dplyr::arrange(pos) ->
  .sort_variant

  .forplot_cluster_cell_variant |>
    dplyr::group_by(cluster, variant) |>
    dplyr::summarise(
      mean_cluster_variant_af = mean(af, na.rm = T),
      sum_cluster_variant_depth = sum(depth, na.rm = T),
      max_cluster_variant_depth = max(depth, na.rm = T)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      sum_cluster_variant_depth_log2 = log2(sum_cluster_variant_depth + 1)
    ) |>
    dplyr::filter(variant %in% .sort_variant$variant) ->
  .cluster_variant_af

  .c <- unique(.cluster_variant_af$cluster)
  .v <- unique(.cluster_variant_af$variant)

  tibble::tibble(
    cluster = rep(.c, each = length(.v)),
    variant = rep(.v, length(.c))
  ) |>
    dplyr::left_join(
      .cluster_variant_af,
      by = c("cluster", "variant")
    ) |>
    dplyr::mutate(
      rect_type = ifelse(max_cluster_variant_depth >= 10, "white", "grey")
    ) |>
    dplyr::mutate(
      variant = factor(variant, .sort_variant$variant)
    ) ->
  .no_depth

  # .theforplot |>
  #   dplyr::group_by(variant) |>
  #   dplyr::summarise(maf = mean(af, na.rm = T)) |>
  #   dplyr::arrange(-maf) ->
  # .sort_variant

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
    dplyr::inner_join(
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
  library(ggnewscale)

  .cl <- as.character(unique(.haplo_forplot$cluster)[[1]])
  .haplo_forplot |>
    dplyr::select(variant, pos, cluster, mean_cluster_variant_af) |>
    dplyr::distinct() |>
    dplyr::filter(cluster == .cl) ->
  .forlabel

  ggplot() +
    ggrepel::geom_text_repel(
      data = .forlabel,
      aes(
        x = pos,
        y = mean_cluster_variant_af,
        label = variant,
      ),
      size = 3,
      # direction = "y",
      # nudge_x = -0.2,
      direction = "y",
      # hjust = "right",
      segment.size = 0,
    ) +
    theme(
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      # axis.line.y.left = element_line(color = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "right",
      legend.key = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), ,
      axis.line.y = element_blank(),
    ) ->
  p_label

  .haplo_forplot |>
    ggplot() +
    ggh4x::facet_wrap2(
      ~cluster,
      ncol = 1,
      strip.position = "right",
      strip = ggh4x::strip_themed(
        background_y = elem_list_rect(
          fill = pcc$color
        ),
        text_y = elem_list_text(
          colour = "white",
          face = c("bold")
        ),
        by_layer_y = FALSE,
      )
    ) +
    # geom_violin(
    #   aes(
    #     x = pos,
    #     y = af,
    #     fill = mean_cluster_variant_af
    #   ),
    #   alpha = 0.5,
    #   size = 1,
    #   color = NA,
    #   show.legend = FALSE
    # ) +
    # scale_fill_gradient2(
    #   name = "AF",
    #   low = "white",
    #   mid = "red",
    #   high = "#3B0049",
    #   midpoint = 0.5,
    # ) +
    ggbeeswarm::geom_quasirandom(
      aes(
        x = pos,
        y = af,
        color = af
      ),
      size = 1,
      dodge.width = .75,
      alpha = .5,
    ) +
    scale_color_gradient2(
      name = "AF",
      low = "white",
      mid = "red",
      high = "#3B0049",
      midpoint = 0.5,
    ) +
    theme(
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line.y.left = element_line(color = "black"),
      # axis.line.x.bottom = element_line(color = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "right",
      legend.key = element_blank(),
      axis.title.y = element_text(color = "black"),
      # axis.title.y = element_blank(),
      axis.text.y = element_text(color = "black"),
      # legend.text = element_text(
      #   size = 14,
      #   color = "black"
      # ),
      # legend.title = element_text(
      #   size = 16,
      #   colour = "black"
      # ),
      # strip.background = element_blank(),
      # strip.text = element_text(
      #   # size = 8,
      #   color = "black",
      #   face = "bold"
      # )
    ) +
    labs(y = "AF") ->
  p_af_cell

  wrap_plots(
    p_label,
    p_af_cell,
    ncol = 1,
    heights = c(1, 8)
  )

  # p_af_cell
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
    file = "heatmap_cell_af.pdf",
    width = 14,
    height = 15
  )
  ComplexHeatmap::draw(object = ch_af_depth$ch_af)
  dev.off()

  pdf(
    file = "heatmap_cell_depth.pdf",
    width = 14,
    height = 15
  )
  ComplexHeatmap::draw(object = ch_af_depth$ch_depth)
  dev.off()
  # log_success("save image")
  log_success("save cell allele heatmap")
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
    file = "heatmap_cluster_af.pdf",
    width = 7,
    height = 15
  )
  ComplexHeatmap::draw(object = cluster_ch_af_depth$ch_af)
  dev.off()

  pdf(
    file = "heatmap_cluster_depth.pdf",
    width = 7,
    height = 15
  )
  ComplexHeatmap::draw(object = cluster_ch_af_depth$ch_depth)
  dev.off()
  log_success("save cluster allele heatmap")
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
log_success("save venn diagram")


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


cmd <- "source {conda_root}/etc/profile.d/conda.sh; conda activate {conda_env}; perl {perlscript} {file.path(jar_path, 'haplogrep3.jar')} {sqlite_path} cell_snvlist.tsv > variant_annotation.tsv" |> glue::glue()
log_debug(cmd)
system(command = cmd)




if (file.exists("variant_annotation.tsv")) {
  cell_anno <- readr::read_tsv("variant_annotation.tsv") |>
    dplyr::mutate(variant = glue::glue("{Position}{Ref}>{Alt}"))

  writexl::write_xlsx(
    x = cell_anno,
    path = "variant_annotation.xlsx"
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

# ! raw --------------------------------------------------------------------

cell_raw_ch_af_depth <- fn_heatmap(
  .forplot = cell_raw_cluster_forplot,
  # .cell_variants = cell_cluster_forplot$forplot$variant,
  .cell_variants = unique(cell_hetero$variant),
  .variant_annotation = variant_annotation
)

{
  pdf(
    file = "heatmap_final_af.pdf",
    width = 25,
    height = 15
  )
  ComplexHeatmap::draw(object = cell_raw_ch_af_depth$ch_af)
  dev.off()

  pdf(
    file = "heatmap_final_depth.pdf",
    width = 14,
    height = 15
  )
  ComplexHeatmap::draw(object = cell_raw_ch_af_depth$ch_depth)
  dev.off()
  log_success("save cluster cell allele heatmap")
}


fn_for_select_af_color <- function() {
  pals <- c("magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo")
  library(grid)
  parallel::mclapply(
    X = pals,
    FUN = function(pal) {
      fn_heatmap(
        .forplot = cell_raw_cluster_forplot,
        .cell_variants = unique(cell_hetero$variant),
        .variant_annotation = variant_annotation,
        col_option = pal,
        show_column_title = TRUE
      ) ->
      .p
      .p$ch_af
    },
    mc.cores = length(pals)
  ) ->
  pl

  pdf(
    file = "af_color_options.pdf",
    width = 28,
    height = 15
  )
  for (p in pl) {
    ComplexHeatmap::draw(object = p)
  }
  dev.off()
}

# fn_for_select_af_color()
# violin plot -------------------------------------------------------------



fn_plot_cell_violin(
  .forplot = cell_raw_cluster_forplot,
  .cell_anno = cell_anno
) -> p_violin

{
  ggsave(
    filename = "violin_final_af.pdf",
    plot = p_violin$p_af,
    device = "pdf",
    width = 24,
    height = 12
  )
  log_success("save cluster cell violin plot")
  ggsave(
    filename = "violin_final_depth.pdf",
    plot = p_violin$p_depth,
    device = "pdf",
    width = 24,
    height = 12
  )
  log_success("save cluster cell violin depth plot")

  writexl::write_xlsx(
    x = list(
      haplo_variant = p_violin$haplo_variant,
      haplo_forplot = p_violin$haplo_forplot
    ),
    path = "violin_tables.xlsx"
  )
  log_success("save cluster cell violin tables")
  data.table::fwrite(p_violin$haplo_variant, "violin_haplo_variant.csv")
  data.table::fwrite(p_violin$haplo_forplot, "violin_haplo_forplot.csv")
  log_success("save cluster cell violin data to CSV files")
}



# ! somatic variant --------------------------------------------------------------------

fn_somatic_variant(
  .haplo_variant = p_violin$haplo_variant,
  .haplo_violin = p_violin$haplo_forplot,
  .n_cells = 10
) -> somatic_variant

readr::write_rds(
  x = somatic_variant,
  file = "variant_somatic.rds"
)
log_success("save somatic variant_somatic.rds")
data.table::fwrite(
  somatic_variant |>
    tibble::enframe() |>
    tidyr::unnest(cols = value) |>
    dplyr::mutate(name = as.character(name)) |>
    dplyr::rename(group = name, value = name),
  "variant_somatic.csv"
)
log_success("save somatic variant_somatic.csv")

parallel::mclapply(
  X = names(somatic_variant),
  FUN = function(.x) {
    .sel_variants <- somatic_variant[[.x]]
    fn_plot_cell_violin(
      .forplot = cell_raw_cluster_forplot,
      .cell_anno = cell_anno,
      .sel_variants = .sel_variants
    ) -> .p_violin

    base_width = 2
    base_height = 4

    width = base_width * ifelse(length(.sel_variants) > 12, 12, length(.sel_variants))
    height = base_height * ceiling(length(.sel_variants) / 12)


    ggsave(
      filename = glue::glue("violin_final_af_{.x}.pdf"),
      plot = .p_violin$p_af,
      device = "pdf",
      width = width,
      height = height
    )

    ggsave(
      filename = glue::glue("violin_final_depth_{.x}.pdf"),
      plot = .p_violin$p_depth,
      device = "pdf",
      width = width,
      height = height
    )
  },
  mc.cores = length(somatic_variant)
)



# ! hotspot --------------------------------------------------------------------

tryCatch(
  expr = {
    p_mtdna <- fn_plot_mtdna()
    p_depth <- fn_plot_coverage(.cluster_coverage = cluster_coverage)

    p_hotspots <- fn_plot_hotspots(
      .forplot = cell_raw_cluster_forplot,
      .cell_anno = cell_anno,
      .sel_variants = somatic_variant$somatic
    )

    ggsave(
      filename = glue::glue("hotspots_final_af_somatic.pdf"),
      plot = wrap_plots(
        p_hotspots,
        p_depth$p_mt_depth_allcell,
        p_mtdna,
        ncol = 1,
        heights = c(1.3, 0.4, 0.1)
      ),
      device = "pdf",
      width = 18,
      height = 9
    )
  },
  error = function(e) {
    log_error(e)
    1
  }
)


# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------

save.image(file = "scMOCHA.rda")
# load(file = "/scr1/users/liuc9/tmp/mito/flu2-a/cromwell-executions/SCMTAH/0138fcd0-c384-42c2-8704-6647767610d2/call-plot_scmtah/execution/scmtah.rda")
