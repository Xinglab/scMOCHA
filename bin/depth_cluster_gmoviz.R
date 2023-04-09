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

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------
fn_plot_coverage <- function(.filename, .celltype) {
  .filename
  .celltype

  bam <- Rsamtools::BamFile(file = .filename)
  Rsamtools::indexBam(bam)

  .coverage <- gmoviz::getCoverage(
    regions_of_interest = "MT",
    bam_file = .filename,
    window_size = 50
  )

  .coverage |>
    data.table::as.data.table() |>
    dplyr::mutate(seqnames = "MT") |>
    plyranges::as_granges() ->
  .coverage_a

  mt_features <- readr::read_rds(
    "/home/liuc9/github/scRNAseq-MitoVariant/fasta/mt_features.grange.gmoviz.rds.gz"
  )


  mt_features |>
    plyranges::filter(type == "Mt_tRNA") ->
  mt_features_pc
  mt_features |>
    plyranges::filter(type != "Mt_tRNA") ->
  mt_features_npc

  mt_ideogram <- gmoviz::getIdeogramData(
    fasta_file = "/home/liuc9/github/scRNAseq-MitoVariant/fasta/rCRS.MT.fasta"
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
    a = purrr::map2(
      .x = filename,
      .y = celltype,
      .f = fn_plot_coverage
    )
  )


# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
