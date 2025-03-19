#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun Apr  9 13:54:35 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------
library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(gmoviz)
# Check if ggtranscript is installed, install if not
if (!requireNamespace("ggtranscript", quietly = TRUE)) {
  message("Installing ggtranscript from GitHub...")
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("dzhang32/ggtranscript")
  library(ggtranscript)
}

# args --------------------------------------------------------------------

args <- commandArgs(TRUE)
mt_features_gmoviz <- args[1]
mt_rcrs_fasta <- args[2]
mt_exons_df <- args[3]
# mt_features_gmoviz <- "/scr1/users/liuc9/mitochondrial/realdata/05-Liming/scmocha-celline/cromwell-executions/scMOCHA/6982bedc-1c08-40aa-8c4f-31b59cebe69b/call-cell_cluster_annotation/inputs/2014965526/mt_features.grange.gmoviz.rds.gz"
# mt_rcrs_fasta <- "/scr1/users/liuc9/mitochondrial/realdata/05-Liming/scmocha-celline/cromwell-executions/scMOCHA/6982bedc-1c08-40aa-8c4f-31b59cebe69b/call-cell_cluster_annotation/inputs/2014965526/rCRS.MT.fasta"
# mt_exons_df <- "/scr1/users/liuc9/mitochondrial/realdata/05-Liming/scmocha-celline/cromwell-executions/scMOCHA/6982bedc-1c08-40aa-8c4f-31b59cebe69b/call-cell_cluster_annotation/inputs/2014965526/mt_exons.df.rds.gz"

# src ---------------------------------------------------------------------

pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
  dplyr::arrange(cancer_types)

# header ------------------------------------------------------------------

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

  mt_features <- readr::read_rds(
    mt_features_gmoviz
  )


  mt_features |>
    plyranges::filter(type == "Mt_tRNA") ->
  mt_features_pc
  mt_features |>
    plyranges::filter(type != "Mt_tRNA") ->
  mt_features_npc

  mt_ideogram <- gmoviz::getIdeogramData(
    fasta_file = mt_rcrs_fasta
  )



  gmoviz::gmovizPlot(
    file_name = glue::glue("gmoviz.{.celltype}.svg"),
    file_type = "svg",
    plotting_functions = {
      gmoviz::gmovizInitialise(
        mt_ideogram,
        space_between_sectors = 25,
        start_degree = 78,
        xaxis_spacing = 30,
        sector_label_size = 1,
        coverage_data = .coverage_a,
        coverage_rectangle = "MT"
      )
      gmoviz::drawFeatureTrack(
        mt_features_npc,
        track_height = 0.13
      )
      gmoviz::drawFeatureTrack(
        mt_features_pc,
        feature_label_cutoff = 80000,
      )
    },
    # legends = legend,
    title = glue::glue("{.celltype} coverage"),
    background_colour = "white",
    width = 10,
    height = 10,
    units = "in"
  )

  .coverage_a
}

# load data ---------------------------------------------------------------
tibble::tibble(
  filename = list.files(
    pattern = "MT_cluster.TAG_CJ_.*.bam$"
  )
) |>
  dplyr::mutate(
    celltype = gsub(
      pattern = "MT_cluster.TAG_CJ_|.bam",
      replacement = "",
      x = filename
    )
  ) ->
celltype_bams


# body --------------------------------------------------------------------

celltype_bams |>
  dplyr::mutate(
    coverage = purrr::map2(
      .x = filename,
      .y = celltype,
      .f = fn_plot_coverage
    )
  ) ->
celltype_bams_cov



celltype_bams_cov |>
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
    celltype,
    pos = end, depth = coverage
  ) ->
coverage


coverage %>%
  ggplot(aes(x = pos, y = depth, fill = celltype)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.03)),
    limits = c(1, 17000),
  ) +
  scale_y_continuous(
    expand = c(0.01, 0),
    label = scales::label_number(scale = 1e-5, suffix = "x10^5")
  ) +
  scale_fill_manual(
    name = "Cell type",
    values = pcc$color,
    guide = guide_legend(nrow = 1)
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    # axis.title.y = element_blank(),
    axis.line.x.bottom = element_line(color = "black"),
    # strip.background = element_rect(fill = NA, colour = "black"),
    strip.background = element_blank(),
    # strip.text = element_text(
    #   color = "black",
    #   face = "bold",
    #   size = 8
    # ),
    strip.text = element_blank(),
    legend.position = "top"
  ) +
  facet_wrap(
    facets = ~celltype,
    ncol = 1,
    strip.position = "top"
  ) +
  labs(y = "Depth") ->
p1


gtf_gene_df <-
  readr::read_rds(
    file = mt_exons_df
  )


gtf_gene_df %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = gene_name
  )) +
  geom_range(aes(fill = transcript_biotype)) +
  geom_intron(
    data = to_intron(gtf_gene_df, "transcript_name"),
    aes(strand = strand)
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.03)),
    limits = c(1, 17000),
    breaks = seq(1000, 17000, 1000),
    labels = seq(1000, 17000, 1000)
  ) +
  # scale_fill_brewer(palette = "Set3")
  ggsci::scale_fill_jama(
    name = "Biotype",
    labels = c("MT rRNA", "MT tRNA", "Protein coding")
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    axis.line = element_line(color = "black"),
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.title.x = element_blank(),
    # axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    x = "Position"
  ) ->
p2

p <- cowplot::plot_grid(
  plotlist = list(p1, p2),
  ncol = 1,
  align = "v",
  rel_heights = c(0.5, 0.5)
)

ggsave(
  filename = "plot-mt-cluster-depth.pdf",
  plot = p,
  device = "pdf",
  width = 17,
  height = 15
)

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(
  file = "depth_cluster_gmoviz.rda"
)
