#!/usr/bin/env Rscript --vanilla
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Nov 12 11:18:42 2024
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
library(prismatic)
library(paletteer)
library(data.table)
# library(rlang)
library(GetoptLong)
library(logger)

# args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.
verbose <- FALSE
spec <- "
Usage: Rscript foorbar.R [options]
Options:

<verbose!> Print messages
"

GetoptLong.options(help_style = "two-column")
GetoptLong(spec, template_control = list(opt_width = 21))

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------
log_threshold(TRACE)
log_layout(layout_glue_colors)

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------
fn_umap_coord <- function(.x) {
  .col_names <- c("UMAP_1", "UMAP_2")

  if ("ref.umap" %in% names(.x@reductions)) {
    .umap <- .x@reductions$ref.umap@cell.embeddings |> data.table::as.data.table()
    colnames(.umap) <- .col_names
    .tsne <- NULL
  } else {
    .umap <- .x@reductions$umap@cell.embeddings |> data.table::as.data.table()
    colnames(.umap) <- .col_names
    .tsne <- .x@reductions$tsne@cell.embeddings |> data.table::as.data.table()
    colnames(.tsne) <- .col_names
  }

  # .umap
  .x@meta.data |>
    dplyr::select(
      celltype
    ) |>
    data.table::as.data.table() ->
  .xx

  .xxx <- dplyr::bind_cols(.umap, .xx) |>
    dplyr::mutate(barcode = rownames(.x@meta.data))

  .xxx
}

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
    dplyr::mutate(
      pos = gsub(
        pattern = ">|[AGCT]",
        "",
        x = variant
      )
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

fn_plot_vaf_featureplot <- function(.thevariant, sc, .cell_annotation = NULL) {
  pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
    dplyr::arrange(cancer_types)
  sc$cell_hetero_coverage |>
    dplyr::filter(variant == .thevariant) ->
  vhc


  sc$umap_coord |>
    dplyr::left_join(vhc, by = "barcode") ->
  vhc_umap

  if (is.null(.cell_annotation)) {
    vhc_umap |>
      ggplot(aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(color = af)) +
      scale_color_gradient2(
        name = "AF",
        low = "grey",
        mid = "gold",
        high = "#F02415"
      ) +
      theme_bw() +
      labs(
        title = .thevariant
      ) +
      theme(
        plot.title = element_text(
          color = "black", face = "bold", hjust = 0.5
        )
      ) ->
    p_feature
    return(p_feature)
  } else {
    vhc_umap |>
      ggplot(aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(color = celltype)) +
      scale_color_manual(
        name = "Cell Type",
        values = pcc$color
      ) +
      theme_bw() +
      labs() +
      theme(
        plot.title = element_blank()
      ) ->
    p_annotation
    return(p_annotation)
  }
}


fn_load_by_path <- function(.dir) {
  library(Seurat)
  azimuth_file <- file.path(.dir, "sc_azimuth.rds.gz")
  cell_hetero_file <- file.path(.dir, "cell.cell_heteroplasmic_df_raw.tsv.gz")
  cell_coverage_file <- file.path(.dir, "cell.coverage.txt.gz")

  sc <- readr::read_rds(azimuth_file)
  sc$umap_coord <- fn_umap_coord(.x = sc$sc_azimuth)
  sc$cell_hetero <- fn_load_hetero(cell_hetero_file)
  sc$cell_coverage <- fn_load_coverage(cell_coverage_file)
  sc$cell_hetero_coverage <- sc$cell_hetero |>
    dplyr::mutate(pos = as.integer(pos)) |>
    dplyr::inner_join(sc$cell_coverage, by = c("barcode", "pos"))

  sc
}

fn_plot_vaf_featureplot_multi <- function(.thevariants, sc) {
  purrr::map(
    .thevariants,
    fn_plot_vaf_featureplot,
    sc = sc
  ) ->
  variant_plots

  fn_plot_vaf_featureplot(.thevariants[[1]], sc, .cell_annotation = TRUE) -> p_annotation

  c(variant_plots, list(p_annotation)) |>
    wrap_plots() +
    guide_area() +
    plot_layout(guides = "collect") &
    theme(
      legend.justification = "left",
      legend.position = "right"
    )
}


# load data ---------------------------------------------------------------


# body --------------------------------------------------------------------
thepath <- "/home/liuc9/github/scMOCHA/06-bigdata/GSE226602/cromwell-executions/scMOCHABatch/192a6bdb-b835-4f39-a21d-9423f9c8165d/call-scMOCHA/shard-13/sub.scMOCHA/c3913f7f-efd1-4d72-9615-2463d684f359/call-gather_outputfiles/execution/GSM7080019"
sc <- fn_load_by_path(thepath)
variant_list_file <- file.path(thepath, "cell_variant_annotation.tsv")
# variant_list <- data.table::fread(variant_list_file)

fn_plot_vaf_featureplot_multi(c("3173G>A", "3176A>T", "3178T>A"), sc) -> p
p
#

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
