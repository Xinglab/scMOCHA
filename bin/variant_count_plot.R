#!/usr/bin/env Rscript --vanilla
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Nov 12 11:35:22 2024
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
    dplyr::mutate(nv = V3 + V4) |>
    dplyr::select(gt = n, pos = V1, group = V2, fw = V3, rv = V4, nv) ->
  cluster_n

  fasta <- Biostrings::readDNAStringSet("/home/liuc9/github/scMOCHA/fasta/rCRS.chrM.fasta")

  fasta$chrM |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "pos") |>
    dplyr::rename(ref = x) |>
    dplyr::mutate(posref = glue::glue("{pos}{ref}")) |>
    dplyr::mutate(pos = as.integer(pos)) ->
  fasta_df

  cluster_n |>
    dplyr::left_join(fasta_df, by = "pos") |>
    dplyr::mutate(gt = factor(gt, levels = c("A", "G", "C", "T"))) |>
    as.data.table() -> cluster_n_temp

  cluster_n_temp[, ratio := nv / sum(nv), by = .(group, pos)]

  cluster_n_temp |>
    dplyr::mutate(
      label = glue::glue("total coverage = {nv} \n forward = {fw}, reverse = {rv} \n ratio = ({round(ratio, 3) * 100}%)")
    ) ->
  cluster_n_forplot

  cluster_n_forplot
}

fn_plot_count <- function(cluster_n_forplot, thepos, group_sel = NA) {
  if (!is.na(group_sel)) {
    cluster_n_forplot |>
      dplyr::filter(group %in% group_sel) ->
    cluster_n_forplot
  }
  if (length(unique(cluster_n_forplot$group)) > 20) {
    stop("The number of unique groups exceeds 20.")
  }

  cluster_n_forplot |>
    dplyr::filter(pos %in% thepos) |>
    dplyr::mutate(pos = as.character(pos)) |>
    ggplot(aes(x = posref, y = gt)) +
    geom_tile(aes(fill = nv)) +
    geom_text(aes(label = label), size = 3.5) +
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
    facet_wrap(~group, ncol = 1, strip.position = "right") ->
  p_tile
  p_tile
}

# load data ---------------------------------------------------------------
thepath <- "/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE157344/cromwell-executions/scMOCHABatch/87f7e7e9-4e27-491a-9125-19a78cddaf64/call-scMOCHA/shard-15/sub.scMOCHA/4d228407-529c-4b74-bc9e-f189ce7b2274/call-gather_outputfiles/execution/GSM4762177"

# body --------------------------------------------------------------------
cluster_n_forplot <- fn_load_count(thepath, type = "cluster")
#
fn_plot_count(cluster_n_forplot, thepos = 217)


# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
