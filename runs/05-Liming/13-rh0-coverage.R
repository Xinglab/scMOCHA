#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed Jul 17 16:56:46 2024
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
fn_plot_coverage <- function(.filename, .celltype) {
  # .filename
  # .celltype
  
  bam <- Rsamtools::BamFile(file = .filename)
  Rsamtools::indexBam(bam)
  
  .coverage <- gmoviz::getCoverage(
    regions_of_interest = "MT",
    bam_file = .filename,
    window_size = 1
  )
  
  .coverage |>
    data.table::as.data.table() |>
    dplyr::mutate(seqnames = "MT") |>
    plyranges::as_granges() ->
    .coverage_a
  
  # mt_features <- readr::read_rds(
  #   mt_features_gmoviz
  # )
  # 
  # 
  # mt_features |>
  #   plyranges::filter(type == "Mt_tRNA") ->
  #   mt_features_pc
  # mt_features |>
  #   plyranges::filter(type != "Mt_tRNA") ->
  #   mt_features_npc
  # 
  # mt_ideogram <- gmoviz::getIdeogramData(
  #   fasta_file = mt_rcrs_fasta
  # )
  
  
  # 
  # gmoviz::gmovizPlot(
  #   file_name = glue::glue("gmoviz.{.celltype}.svg"),
  #   file_type = "svg",
  #   plotting_functions = {
  #     gmoviz::gmovizInitialise(
  #       mt_ideogram,
  #       space_between_sectors = 25,
  #       start_degree = 78,
  #       xaxis_spacing = 30,
  #       sector_label_size = 1,
  #       coverage_data = .coverage_a,
  #       coverage_rectangle = "MT"
  #     )
  #     gmoviz::drawFeatureTrack(
  #       mt_features_npc,
  #       track_height = 0.13
  #     )
  #     gmoviz::drawFeatureTrack(
  #       mt_features_pc,
  #       feature_label_cutoff = 80000,
  #     )
  #   },
  #   # legends = legend,
  #   title = glue::glue("{.celltype} coverage"),
  #   background_colour = "white",
  #   width = 10,
  #   height = 10,
  #   units = "in"
  # )
  
  .coverage_a
}

fn_plot_gene <- function() {
  mt_exons_df <- "/scr1/users/liuc9/mitochondrial/realdata/05-Liming/scmocha-mixed-cellline-high-depth/cromwell-executions/scMOCHA/023d7328-9097-4e50-8c11-19f860c5519e/call-cellranger_count/inputs/2014965526/mt_exons.df.rds.gz"
  
  
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
    theme_genes() +
    theme(
      legend.position = "bottom",
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 14),
      legend.text = element_text(size = 14)
    ) ->
    pg;pg
}

fn_plot_gggene <- function(.coverage) {
  
  pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
    dplyr::arrange(cancer_types)
  
  .coverage |> 
    dplyr::select(pos = pos, group = celltype, nv = depth) |> 
    dplyr::filter(group != "cluster_4") |> 
    dplyr::mutate(
      group2 = plyr::revalue(x = group, replace = c("cluster_0" = "A549", "cluster_1" = "WAL2A", "cluster_2" = "143B", "cluster_3" = "HEK293"))
    ) |> 
    ggplot(aes(x = pos, y = nv)) +
    geom_line(aes(color = group2 )) +
    # geom_vline(xintercept = thepos, color = "red") +
    scale_x_continuous(
      expand = expansion(mult = c(0.01, 0)),
      limits = c(1, 17000),
      breaks = seq(0, 17000, 2000),
      labels = seq(0, 17000, 2000)
    ) +
    scale_y_continuous(
      expand = c(0.01, 0),
      label = scales::label_number()
    ) +
    # ggsci::scale_color_jco(
    scale_color_brewer(
      name = "Dataset",
      palette = "Set1"
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
      axis.title.y = element_text(size = 16, color = "black"),
      axis.text.y = element_text(size = 14, color = "black"),
      legend.text = element_text(
        size = 14,
        color = "black"
      ),
      legend.title = element_text(
        size = 16,
        colour = "black"
      )
    ) +
    labs(y = "Depth") ->
    pp;pp
  
  # wrap_plots(
  #   pp,
  #   pg,
  #   ncol = 1,
  #   heights = c(0.9, 0.1)
  # ) ->
  #   pg_merged_p_read_depth;pg_merged_p_read_depth
}
# load data ---------------------------------------------------------------

# body --------------------------------------------------------------------



# Rh0 ---------------------------------------------------------------------
rh0_path <- "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth/cromwell-executions/scMOCHA/b86127e0-3759-4d98-96b1-22357e8b6a21/call-cell_cluster_annotation/execution"

tibble::tibble(
  filename = list.files(
    path = rh0_path,
    pattern = "MT_cluster.TAG_CJ_.*.bam$"
  )
) |>
  dplyr::mutate(
    celltype = gsub(
      pattern = "MT_cluster.TAG_CJ_|.bam",
      replacement = "",
      x = filename
    )
  ) |> 
  dplyr::mutate(
    coverage = purrr::map2(
      .x = filename,
      .y = celltype,
      .f = \(.filename, .celltype) {
        .path <- file.path(rh0_path, .filename)
        fn_plot_coverage(.path, .celltype)
      }
    )
  ) |> 
  dplyr::mutate(
    coverage = purrr::map(
      .x = coverage,
      .f = data.table::as.data.table
    )
  ) |>
  dplyr::select(-filename) |>
  dplyr::mutate(celltype = factor(celltype)) |>
  tidyr::unnest(cols = coverage) |>
  dplyr::select(
    celltype, pos = end, depth = coverage
  ) ->
  coverage_rh0




p_rh0

# WT ----------------------------------------------------------------------
wt_path <- "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth/cromwell-executions/scMOCHA/0f65a458-3982-4c90-b85e-83dd074a99ce/call-cell_cluster_annotation/execution"

tibble::tibble(
  filename = list.files(
    path = wt_path,
    pattern = "MT_cluster.TAG_CJ_.*.bam$"
  )
) |>
  dplyr::mutate(
    celltype = gsub(
      pattern = "MT_cluster.TAG_CJ_|.bam",
      replacement = "",
      x = filename
    )
  ) |> 
  dplyr::mutate(
    coverage = purrr::map2(
      .x = filename,
      .y = celltype,
      .f = \(.filename, .celltype) {
        .path <- file.path(wt_path, .filename)
        fn_plot_coverage(.path, .celltype)
      }
    )
  ) |> 
  dplyr::mutate(
    coverage = purrr::map(
      .x = coverage,
      .f = data.table::as.data.table
    )
  ) |>
  dplyr::select(-filename) |>
  dplyr::mutate(celltype = factor(celltype)) |>
  tidyr::unnest(cols = coverage) |>
  dplyr::select(
    celltype, pos = end, depth = coverage
  ) ->
  coverage_wt


p_ggene <- fn_plot_gene()
fn_plot_gggene(coverage) ->p_rh0
fn_plot_gggene(coverage_wt) ->p_wt

p_wt / p_rh0 / p_ggene + plot_layout(
  heights = c(8,8,1)
) ->
  p_merge;p_merge
ggplot2::ggsave(
  filename = "combined_read_depth_WT_Rh0.pdf",
  plot = p_merge,
  path = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth",
  width = 15,
  height = 8
)
# 
# p_wt 
# p_rh0
# ggplot2::ggsave(
#   filename = "combined_read_depth_WT.pdf",
#   plot = p_wt,
#   path = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth",
#   width = 15,
#   height = 7
# )
# ggplot2::ggsave(
#   filename = "combined_read_depth_Rh0.pdf",
#   plot = p_rh0,
#   path = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth",
#   width = 15,
#   height = 7
# )



# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------