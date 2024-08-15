#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Mon Jul 15 14:50:21 2024
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
library(prismatic)
library(paletteer)
library(data.table)
#library(rlang)
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

# log_info('Starting the script...')
# log_debug('This is the second log line')
# log_trace('Note that the 2nd line is being placed right after the 1st one.')
# log_success('Doing pretty well so far!')
# log_warn('But beware, as some errors might come :/')
# log_error('This is a problem')
# log_debug('Note that getting an error is usually bad')
# log_error('This is another problem')
# log_fatal('The last problem')

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

fn_plot_vaf_featureplot<- function(.thevariant, cell_hetero_coverage, umap_coord) {
  cell_hetero_coverage |> 
    dplyr::filter(variant == .thevariant) ->
    vhc
  
  
  umap_coord|> 
    dplyr::left_join(vhc, by = "barcode") ->
    vhc_umap
  
  vhc_umap |> 
    ggplot(aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = af)) +
    scale_color_gradient2(
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
        color = "black",  face = "bold", hjust = 0.5
      )
    )
}

# load data ---------------------------------------------------------------
library(Seurat)
azimuth_file <- "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth2/cromwell-executions/scMOCHA/139358d8-df39-4274-b931-9c42b8d9c3bb/call-gather_outputfiles/execution/WT/sc_azimuth.rds.gz"
cell_hetero_file <- "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth2/cromwell-executions/scMOCHA/139358d8-df39-4274-b931-9c42b8d9c3bb/call-gather_outputfiles/execution/WT/cell.cell_heteroplasmic_df_raw.tsv.gz"
cell_coverage_file <- "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth2/cromwell-executions/scMOCHA/139358d8-df39-4274-b931-9c42b8d9c3bb/call-gather_outputfiles/execution/WT/cell.coverage.txt.gz"



# body --------------------------------------------------------------------
sc <- readr::read_rds(azimuth_file)
# hetero <- data.table::fread(hetero_file)


sc$umap_coord <- fn_umap_coord(.x = sc$sc_azimuth)
cell_hetero <- fn_load_hetero(cell_hetero_file)
cell_coverage <- fn_load_coverage(cell_coverage_file)
cell_hetero_coverage <- cell_hetero |> 
  dplyr::mutate(pos = as.integer(pos)) |> 
  dplyr::inner_join(cell_coverage, by = c("barcode", "pos"))



 # 143B
variant_list_file <- "/mnt/isilon/u01_project/PT/Comparison_from_different_cutoff_combinations_0708/Comparison_result_500_8000_2nd_C0.05_S0.2/Chunjie_specific_143B.txt"
variant_list <- data.table::fread(variant_list_file)

variant_list |> 
  dplyr::mutate(
    p = purrr::map(
      .x = Variant, .f = fn_plot_vaf_featureplot,
      cell_hetero_coverage, umap_coord = sc$umap_coord
    )
  ) ->
  variant_list_p


wrap_plots(
  variant_list_p$p
) +
  guide_area() +
  plot_layout(guides = "collect") ->p;p

ggsave(
  filename = "143B-featureplot.pdf",
  path = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth/featureplot",
  plot = p,
  width = 9, 
  height = 5
  
)


tibble::tibble(
  Variant = "6776T>C"
) |> 
  dplyr::mutate(
    p = purrr::map(
      .x = Variant, .f = fn_plot_vaf_featureplot,
      cell_hetero_coverage, umap_coord = sc$umap_coord
    )
  ) ->
  p;p

p$p


vhc_umap |> 
  dplyr::filter(!is.na(variant)) ->
  vhc_umap_filter


vhc_umap_filter |> 
  dplyr::arrange(celltype, -af)  ->
  thesort

vhc_umap_filter |> 
  dplyr::mutate(barcode = factor(barcode, thesort$barcode)) |> 
  ggplot(aes(
    x = barcode,
    y = variant
  )) +
  geom_tile(
    aes(fill = af)
  ) +
  scale_fill_gradient(
    low = "white",
    high = "red"
  ) +
  scale_y_discrete(
    expand = c(0, 0)
  ) +
  theme(
    panel.background = element_rect(
      color = "black",
      fill = NA
    ),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 16, color = 'black'),
    axis.title.y = element_text(size = 18, color = "black")
  ) +
  labs(
    y = "AF"
  ) ->
  v_p_af;v_p_af

vhc_umap_filter |> 
  dplyr::mutate(barcode = factor(barcode, thesort$barcode)) |> 
  ggplot(aes(
    x = barcode,
    y = variant
  )) +
  geom_tile(
    aes(fill = depth)
  ) +
  scale_fill_gradient(
    low = "#FDE725FF",
    high = "#440154FF"
  ) +
  scale_y_discrete(
    expand = c(0, 0)
  ) +
  theme(
    panel.background = element_rect(
      color = "black",
      fill = NA
    ),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 16, color = 'black'),
    axis.title.y = element_text(size = 18, color = "black")
  ) +
  labs(
    y = "DP"
  ) ->
  v_p_dp;v_p_dp

vhc_umap_filter |> 
  dplyr::mutate(barcode = factor(barcode, thesort$barcode)) |>
  ggplot(aes(
    x = barcode, 
    y = variant
  )) +
  geom_tile(
    aes(fill = celltype)
  ) +
  scale_fill_brewer(palette = "Set3") +
  scale_y_discrete(
    expand = c(0, 0)
  ) +
  theme(
    panel.background = element_rect(
      color = "black",
      fill = NA
    ),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 16, color = 'black'),
    axis.title.y = element_text(size = 18, color = "black")
  ) +
  labs(
    y = "Cell type"
  ) ->
  v_p_ct;v_p_ct

wrap_plots(
  v_p_dp, v_p_af , v_p_ct
) +
plot_layout(
  heights = c(10, 10, 2),
  guides = "collect"
) ->
  p_collect;p_collect
ggsave(
  filename = "variant_6776_cell_level.pdf",
  plot = p_collect,
  path = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth2",
  height = 7,
  width = 10
  
)


# cell level 

tibble::tibble(
  path = list.files(
    "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth2/cromwell-executions/scMOCHA/139358d8-df39-4274-b931-9c42b8d9c3bb/call-call_mt_variants/execution/cluster/final",
    "*cluster.*.txt.gz*",
    full.names = T
  )
) |> 
  dplyr::mutate(d = purrr::map(path, data.table::fread)) |> 
  dplyr::mutate(n = basename(path)) |> 
  dplyr::mutate(n = gsub("cluster.|.txt.gz", "", n)) |> 
  dplyr::select(n, d) ->
  cluster

cluster |> 
  dplyr::filter(n != "coverage") |> 
  tidyr::unnest(cols = d) |> 
  dplyr::mutate(nv = V3 + V4) |> 
  dplyr::select(gt = n, pos = V1, group = V2, nv) ->
  cluster_n

cluster |> 
  dplyr::filter(n != "coverage") |> 
  tidyr::unnest(cols = d) |> 
  dplyr::mutate(nv = V3 + V4) |> 
  dplyr::select(gt = n, pos = V1, group = V2, fw = V3, rv = V4, nv) ->cluster_n

fasta <- Biostrings::readDNAStringSet("/home/liuc9/github/scMOCHA/fasta/rCRS.chrM.fasta")

fasta$chrM |> as.data.frame() |> 
  tibble::rownames_to_column(var = "pos") |> 
  dplyr::rename(ref = x) |> 
  dplyr::mutate(posref = glue::glue("{pos}{ref}")) |> 
  dplyr::mutate(pos = as.integer(pos)) ->
  fasta_df

cluster_n |> 
  dplyr::left_join(fasta_df, by = "pos") |> 
  # dplyr::mutate(pos = as.character(pos)) |> 
  dplyr::mutate(gt = factor(gt, levels = c("A", "G", "C", "T"))) |> 
  dplyr::group_by(group, pos) |> 
  # dplyr::group_by(pos, gt) |> 
  dplyr::mutate(ratio = nv / sum(nv)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(
    label = glue::glue("total coverage = {nv} \n forward = {fw}, reverse = {rv} \n ratio = ({round(ratio, 3) * 100}%)")
  ) |> 
  dplyr::mutate(
    group2 = plyr::revalue(x = group, replace = c("cluster_0" = "WAL2A-1", "cluster_1" = "WAL2A-2", "cluster_2" = "HEK293", "cluster_3" = "A549", "cluster_4" = "143B"))
  ) ->
  cluster_n_forplot
thepos = 6776

cluster_n_forplot |> 
  dplyr::filter(pos %in% thepos) |> 
  dplyr::mutate(pos = as.character(pos)) |> 
  ggplot(aes(x = posref, y = gt)) +
  geom_tile(aes(fill = nv)) +
  geom_text(aes(label= label), size = 3.5) +
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
  facet_wrap(~group2, ncol = 1, strip.position = "right") ->
  p_tile;p_tile

# A549
variant_list_file <- "/mnt/isilon/u01_project/PT/Comparison_from_different_cutoff_combinations_0708/Comparison_result_500_8000_2nd_C0.05_S0.2/Chunjie_specific_A549.txt"
variant_list <- data.table::fread(variant_list_file)

variant_list |> 
  dplyr::mutate(
    p = purrr::map(
      .x = Variant, .f = fn_plot_vaf_featureplot,
      cell_hetero_coverage, umap_coord = sc$umap_coord
    )
  ) ->
  variant_list_p



wrap_plots(
  variant_list_p$p
) +
  guide_area() +
  plot_layout(guides = "collect") ->p;p
ggsave(
  filename = "A549-featureplot.pdf",
  path = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth/featureplot",
  plot = p,
  width = 12, 
  height = 7
  
)


# HEK293
variant_list_file <- "/mnt/isilon/u01_project/PT/Comparison_from_different_cutoff_combinations_0708/Comparison_result_500_8000_2nd_C0.05_S0.2/Chunjie_specific_HEK293.txt"
variant_list <- data.table::fread(variant_list_file)

variant_list |> 
  dplyr::mutate(
    p = purrr::map(
      .x = Variant, .f = fn_plot_vaf_featureplot,
      cell_hetero_coverage, umap_coord = sc$umap_coord
    )
  ) ->
  variant_list_p



wrap_plots(
  variant_list_p$p,
  ncol = 4
) +
  guide_area() +
  plot_layout(guides = "collect") ->p;p
ggsave(
  filename = "HEK293-featureplot.pdf",
  path = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth/featureplot/",
  plot = p,
  width = 12, 
  height = 13
  
)



# WAL2A
variant_list_file <- "/mnt/isilon/u01_project/PT/Comparison_from_different_cutoff_combinations_0708/Comparison_result_500_8000_2nd_C0.05_S0.2/Chunjie_specific_WAL2A.txt"
variant_list <- data.table::fread(variant_list_file)

variant_list |> 
  dplyr::mutate(
    p = purrr::map(
      .x = Variant, .f = fn_plot_vaf_featureplot,
      cell_hetero_coverage, umap_coord = sc$umap_coord
    )
  ) ->
  variant_list_p



wrap_plots(
  variant_list_p$p,
  ncol = 3
) +
  guide_area() +
  plot_layout(guides = "collect") ->p;p
ggsave(
  filename = "WAL2A-featureplot.pdf",
  path = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth/featureplot/",
  plot = p,
  width = 9, 
  height = 5
  
)

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------