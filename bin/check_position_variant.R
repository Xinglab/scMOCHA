#!/usr/bin/env Rscript --vanilla
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: `r date()`
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

# future: :plan(future: :multisession, workers = 10)

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

  data.table::fread(input = .filename) -> .d

  data.table::setnames(.d, "V1", "barcode")
  .d <- data.table::melt(.d, id.vars = "barcode", variable.name = "variant", value.name = "af")
  .d[, pos := gsub(pattern = ">|[AGCT]", replacement = "", x = variant)]
  .d[, pos := as.integer(pos)]
  .d
}

fn_load_coverage <- function(.filename) {
  data.table::fread(
    input = .filename,
    sep = ",",
    col.names = c("pos", "barcode", "depth")
  ) ->
  .d
  .d[, depth := log2(depth + 1)]
  .d
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


fn_load_by_path <- function(thepath) {
  library(Seurat)
  azimuth_file <- file.path(thepath, "sc_azimuth.rds.gz")
  cell_hetero_file <- file.path(thepath, "cell.cell_heteroplasmic_df_raw.tsv.gz")
  cell_coverage_file <- file.path(thepath, "cell.coverage.txt.gz")

  sc <- readr::read_rds(azimuth_file)
  sc$umap_coord <- fn_umap_coord(.x = sc$sc_azimuth)
  sc$cell_hetero <- fn_load_hetero(cell_hetero_file)
  sc$cell_coverage <- fn_load_coverage(cell_coverage_file)
  sc$cell_hetero_coverage <- sc$cell_hetero |>
    dplyr::left_join(sc$cell_coverage, by = c("barcode", "pos"))

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
    wrap_plots(ncol = 4) +
    guide_area() +
    plot_layout(guides = "collect") &
    theme(
      legend.justification = "left",
      legend.position = "right"
    )
}

fn_load_count <- function(thepath, type = c("cluster", "cell")) {
  type <- match.arg(type)

  pattern <- if (type == "cluster") {
    "*cluster.*.txt.gz*"
  } else {
    "*cell.*.txt.gz*"
  }

  tibble::tibble(
    path = list.files(
      thepath,
      pattern,
      full.names = T
    )
  ) |>
    dplyr::filter(!grepl("coverage", x = path)) |>
    dplyr::mutate(d = purrr::map(path, data.table::fread)) |>
    dplyr::mutate(n = basename(path)) |>
    dplyr::mutate(n = gsub(paste0(type, ".|.txt.gz"), "", n)) |>
    dplyr::select(n, d) |>
    tidyr::unnest(cols = d) |>
    as.data.table() |>
    dplyr::mutate(nv = V3 + V4) |>
    dplyr::select(gt = n, pos = V1, group = V2, fw = V3, rv = V4, nv) ->
  cluster_n

  fasta <- Biostrings::readDNAStringSet("/home/liuc9/github/scMOCHA/fasta/rCRS.chrM.fasta")

  fasta$chrM |>
    as.data.table() |>
    tibble::rownames_to_column(var = "pos") |>
    dplyr::rename(ref = x) |>
    dplyr::mutate(posref = glue::glue("{pos}{ref}")) |>
    dplyr::mutate(pos = as.integer(pos)) ->
  fasta_df

  cluster_n |>
    dtplyr::lazy_dt() |>
    dplyr::left_join(fasta_df, by = "pos") |>
    dplyr::mutate(gt = factor(gt, levels = c("A", "G", "C", "T"))) |>
    as.data.table() ->
  cluster_n_temp

  cluster_n_temp[, ratio := nv / sum(nv), by = .(group, pos)]

  cluster_n_temp |>
    dplyr::mutate(
      label = glue::glue("total coverage = {nv} \n forward = {fw}, reverse = {rv} \n ratio = ({round(ratio, 3) * 100}%)")
    ) ->
  cluster_n_forplot

  cluster_n_forplot
}

fn_plot_count <- function(cluster_n_forplot, thepos, group_sel = NA) {
  pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
    dplyr::arrange(cancer_types)
  if (!all(is.na(group_sel))) {
    cluster_n_forplot |>
      dplyr::filter(group %in% group_sel) ->
    cluster_n_forplot
  }
  if (length(unique(cluster_n_forplot$group)) > 10) {
    # stop("The number of unique groups exceeds 50.")
    cluster_n_forplot |>
      dplyr::filter(pos == thepos) |>
      dplyr::arrange(-nv) |>
      dplyr::slice(1:10) ->
    cluster_n_forplot
  }

  cluster_n_forplot |>
    dplyr::filter(pos %in% thepos) |>
    dplyr::mutate(pos = as.character(pos)) |>
    ggplot(aes(x = posref, y = gt)) +
    geom_tile(aes(fill = nv)) +
    geom_text(aes(label = label)) +
    scale_fill_gradient(
      low = "white",
      high = "red"
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      # axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(
        color = "black",
        size = 18
      ),
      legend.position = "none ",
      plot.title = element_text(
        size = 16,
        hjust = 0.5
      ),
      strip.background = element_rect(
        fill = NA,
        color = "black",
      ),
      strip.text = element_text(
        color = "black",
        size = 14,
        face = "bold"
      ),
      axis.line = element_line(
        color = "black"
      )
    ) +
    ggh4x::facet_wrap2(
      ~group,
      ncol = 8,
      strip.position = "top",
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = pcc$color
        ),
        text_x = ggh4x::elem_list_text(
          colour = "white",
          face = c("bold")
        ),
        by_layer_y = FALSE,
      )
    ) ->
  p_tile
  p_tile
}

fn_plot_count_multi <- function(cluster_n_forplot, thepos, group_sel = NA) {
  theposes |>
    purrr::map(~ fn_plot_count(cluster_n_forplot, thepos = ., group_sel = NA)) |>
    wrap_plots(ncol = 1) +
    guide_area() +
    plot_layout(guides = "collect") &
    theme(
      legend.justification = "left",
      legend.position = "none"
    ) ->
  p_count
  p_count
}

# load data ---------------------------------------------------------------

thepath <- "/home/liuc9/github/scMOCHA/06-bigdata/GSE226602/cromwell-executions/scMOCHABatch/192a6bdb-b835-4f39-a21d-9423f9c8165d/call-scMOCHA/shard-13/sub.scMOCHA/c3913f7f-efd1-4d72-9615-2463d684f359/call-gather_outputfiles/execution/GSM7080019"

# variant_list_file <- file.path(thepath, "cell_variant_annotation.tsv")

# load sc
sc <- fn_load_by_path(thepath)
# load count
cluster_n_forplot <- fn_load_count(thepath, type = "cluster")


# targeted variant and position ------------------------------------------
thevariants <- c(
  "2191A>C", "2192A>T", "2193T>A",
  "3173G>A", "3176A>T", "3178T>A"
)
theposes <- thevariants |>
  purrr::map(~ gsub(pattern = "[>|AGCT]", "", x = .)) |>
  purrr::map_int(as.integer)

fn_plot_vaf_featureplot_multi(
  .thevariants = thevariants,
  sc = sc
) -> p_vaf_feature
p_vaf_feature

fn_plot_count_multi(
  cluster_n_forplot,
  thepos = theposes
) -> p_count
p_count

# body --------------------------------------------------------------------


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
