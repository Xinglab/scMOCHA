# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu May 25 23:37:36 2023
# @DESCRIPTION:

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
metafile <- "/scr1/users/liuc9/mitochondrial/realdata/04-HubMAP/hubmap-datasets-metadata-2023-05-26_03-32-34.tsv"

# body --------------------------------------------------------------------

metadata <- readr::read_tsv(file = metafile) |>
  dplyr::filter(uuid != "#")

metadata |>
  dplyr::glimpse()

metadata |>
  dplyr::count(assay_type) |>
  dplyr::arrange(assay_type) |>
  print(n = Inf)

metadata |>
  dplyr::filter(
    grepl(
      pattern = "scrnaseq",
      x = assay_type,
      ignore.case = TRUE
    )
  ) ->
  metadata_scrnaseq

metadata_scrnaseq |>
  dplyr::glimpse()
metadata_scrnaseq |>
  dplyr::select(
    assay_type,
    mapped_consortium,
    origin_sample.mapped_organ
  ) ->
  metadata_scrnaseq_sel


metadata_scrnaseq_sel |>
  dplyr::mutate(
    assay_type = ifelse(
      assay_type == "scRNAseq-10xGenomics",
      "scRNAseq-10xGenomics-v3",
      assay_type
    )
  ) |>
  ggplot(aes(
    x = forcats::fct_infreq(origin_sample.mapped_organ),
    # x = origin_sample.mapped_organ,
    fill = assay_type
  )) +
  geom_bar(
    # width = 0.7
  ) +
  geom_text(stat='count', aes(label=..count..), vjust=0,size = 10) +
  ggsci::scale_fill_jama(
    name = "Seq"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1), add = 0)
  ) +
  scale_x_discrete(
    expand = expansion(mult = 0.1, add = 0)
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(
      color = "black",
      size = 18,
      face = "bold"
    ),
    axis.text = element_text(
      color = "black",
      size = 18,
      face = "bold"
    ),
    axis.title.x = element_blank(),
    legend.position = c(0.8, 0.9)
  ) +
  labs(
    y = "Count"
  ) ->
  pbar

ggsave(
  filename = "hubmap-scrnaseq-data.pdf",
  plot = pbar,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA/04-HubMAP/output",
  width = 11,
  height = 6
)

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------
