#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed Mar 22 17:19:46 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(ggtranscript)
library(GetoptLong)
library(logger)

# src ---------------------------------------------------------------------
# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.

# depthfile <- "WT/outs/possorted_genome_bam.MT.depth"
# outfile <- "WT/outs/possorted_genome_bam.MT.depth.pdf"
# mt_exons_df <- "/scr1/users/liuc9/mitochondrial/realdata/05-Liming/scmocha-mixed-cellline-high-depth/cromwell-executions/scMOCHA/023d7328-9097-4e50-8c11-19f860c5519e/call-cellranger_count/inputs/2014965526/mt_exons.df.rds.gz"

verbose <- FALSE

spec <- "
Usage: Rscript scMOCHA.R [options]

Options:
<depthfile=s> possorted_genome_bam.MT.depth
<outfile=s> possorted_genome_bam.MT.depth.pdf
<mt_exons_df=s> mt_exons.df.rds.gz
<verbose!> Print messages
"

GetoptLong.options(help_style = "two-column")
GetoptLong(spec, template_control = list(opt_width = 50))

# args <- commandArgs(TRUE)
# 
# depthfile <- args[1]
# outfile <- args[2]
# mt_exons_df <- args[3]

# depthfile <- "/scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/flu2/Flu2/outs/MT.depth"
# mt_exons_df <- "/home/liuc9/github/scMOCHA/fasta/mt_exons.df.rds.gz"


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
coverage <- data.table::fread(
  input = depthfile,
  col.names =  c("chr", "pos", "depth")
)

# body --------------------------------------------------------------------
# conn <- DBI::dbConnect(
#   duckdb::duckdb(),
#   "/mnt/isilon/xing_lab/liuc9/refdata/ensembl/Homo_sapiens.GRCh38.107.gtf.plyranges.duckdb"
# )
#
# gtf_gene <- dplyr::tbl(conn, "grch38_107_plyranges") |>
#   dplyr::filter(seqnames == "MT") |>
#   data.table::as.data.table()
#
# DBI::dbDisconnect(conn,  shutdown=TRUE)
#
#
# gtf_gene |>
#   dplyr::filter(type == "exon") |>
#   as.data.frame() ->
#   gtf_gene_df
# readr::write_rds(
#   x = gtf_gene_df,file = "/home/liuc9/github/scMOCHA/fasta/mt_exons.df.rds.gz"
# )

gtf_gene_df <-
  readr::read_rds(
    file = mt_exons_df
  )

coverage %>%
  ggplot(aes(x=pos, y = depth)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.03)),
    limits = c(1, 17000),
  ) +
  scale_y_continuous(
    expand = c(0.01, 0)
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  labs(y = "Depth") ->
  p1


gtf_gene_df %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = gene_name
  )) +
  geom_range( aes(fill = transcript_biotype)) +
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
  rel_heights = c(0.3, 0.7)
)

ggsave(
  filename = outfile,
  plot = p,
  device = "pdf",
  width = 10,
  height = 13
)

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------