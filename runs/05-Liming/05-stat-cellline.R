#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed Dec 13 15:24:38 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
library(data.table)
#library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
datadir <- "/home/liuc9/github/scMOCHA/05-Liming/cellline/torun"

# body --------------------------------------------------------------------

tibble::tibble(
  torun = list.dirs(
    path = datadir,
    full.names = T,
    recursive = F 
  )
) |> 
  dplyr::mutate(
    log = purrr::map_chr(
      .x = torun,
      .f = \(.x) {
        .log <- "{basename(.x)}.log" |> glue::glue()
        file.path(
          .x,
          .log
        )
      }
    )
  ) |> 
  dplyr::mutate(
    projectname = basename(torun)
  ) ->
  logfile

logfile |> 
  dplyr::mutate(
    outputdir = purrr::map_chr(
      .x = log,
      .f = \(.log) {
        # .log <- logfile$log[[1]]
        
        .l <- readr::read_lines(
          file = .log
        ) 
        
        .l |> 
          stringr::str_detect(
            "scMOCHA.output_dir_tar_gz"
          ) |> 
          which() ->
          .arr
        
        if(length(.arr) == 0) {return(NA_character_)}
        .l[[.arr[[1]]]] |> 
          gsub("\"| |,|.tar.gz", "", x = _) |> 
          strsplit(
            split = ":"
          ) ->
          .s
        .s[[1]][[2]]
      }
    )
  ) |> 
  dplyr::select(
    projectname, outputdir
  ) |> 
  dplyr::filter(!is.na(outputdir)) ->
  outdir

readr::write_rds(
  x = outdir,
  file = "/home/liuc9/github/scMOCHA/05-Liming/cellline/torun/outdir.rds.gz"
)

future::plan(future::multisession, workers = 7)
outdir |> 
  dplyr::mutate(
    cluster = furrr::future_map(
      .x = outputdir,
      .f = \(.x) {
        # .x <- outdir$outputdir[[1]]
        qc_cell_stats <- readxl::read_xlsx(
          path = file.path(
            .x, 
            "qc_cell_stats.xlsx"
          )
        )
        read_depth <- data.table::fread(
          input = file.path(
            .x,
            "possorted_genome_bam.MT.depth"
          )
        )
        celltype_ratio <- readr::read_tsv(
          file.path(
            .x,
            "celltype_ratio.tsv"
          )
        )
        # cell_variant_annotation <- readr::read_tsv(
        #   file.path(
        #     .x,
        #     "cell_variant_annotation.tsv"
        #   )
        # )
        # 
        # cell_heteroplasmic_df_raw <- readr::read_tsv(
        #   file.path(
        #     .x,
        #     "cell.cell_heteroplasmic_df_raw.tsv.gz"
        #   )
        # )
        # cell_coverage <- readr::read_tsv(
        #   file.path(
        #     .x,
        #     "cell.coverage.txt.gz"
        #   )
        # )
        
        tibble::tibble(
          qc_cell_stats = list(qc_cell_stats),
          read_depth = list(read_depth),
          celltype_ratio = list(celltype_ratio),
        #   cell_variant_annotation = list(cell_variant_annotation),
        #   cell_heteroplasmic_df_raw = list(cell_heteroplasmic_df_raw),
        #   cell_coverage = list(cell_coverage)
        )
      }
    )
  )  ->
  alldataloaded
future::plan(future::sequential)

# qc_cell_stats -----------------------------------------------------------


alldataloaded |> 
  tidyr::unnest(cols = cluster) |> 
  dplyr::select(
    projectname,
    qc_cell_stats
  ) |> 
  tidyr::unnest(cols = qc_cell_stats) ->
  qc_cell_stats
writexl::write_xlsx(
  x = qc_cell_stats,
  path = "/home/liuc9/github/scMOCHA/05-Liming/cellline/torun/qc_cell_stats.xlsx"
)


alldataloaded |> 
  tidyr::unnest(cols = cluster) |> 
  dplyr::select(projectname, celltype_ratio) |> 
  tidyr::unnest(cols = celltype_ratio) |> 
  ggplot(aes(
    x = projectname,
    y = ratio
  )) +
  geom_col(aes(fill = celltype)) +
  ggsci::scale_fill_lancet(name = "Clusters") +
  theme_bw() +
  labs(
    x = "",
    y = "Ratio"
  ) -> 
  p_ratio;p_ratio

ggsave(
  filename = "cluster_ratio.pdf",
  plot = p_ratio,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA/05-Liming/cellline/torun",
  width = 8,
  height = 5
)

# read depth --------------------------------------------------------------

rh0_depth <- readr::read_tsv("/home/liuc9/github/scMOCHA/05-Liming/scmocha-celline/cromwell-executions/scMOCHA/4adb4525-bf49-445d-b433-40c4b84a9932/call-cellranger_count/execution/GEX_Rh0/outs/possorted_genome_bam.MT.depth")

alldataloaded |> 
  tidyr::unnest(cols = cluster) |> 
  dplyr::glimpse()


alldataloaded |> 
  tidyr::unnest(cols = cluster) |> 
  dplyr::select(projectname, read_depth) |> 
  tidyr::unnest(cols = read_depth) ->
  read_depth

read_depth |> 
  dplyr::mutate(
    V3 = V3 / 1000000
  ) |> 
  ggplot(aes(
    x = V2,
    y = V3
  )) +
  geom_line(
    aes(color = projectname),
    stat = "identity"
  ) +
  geom_vline(xintercept = 3243, color = "black") +
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0)),
    limits = c(1, 17000),
    breaks = seq(0, 17000, 2000),
    labels = seq(0, 17000, 2000)
  ) +
  scale_y_continuous(
    expand = c(0.01, 0)
  ) +
  ggsci::scale_color_aaas(
    name = "Dataset",
    guide = guide_legend(nrow = 1)
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.key = element_blank()
  ) +
  labs(y = "Depth (x10^6)") ->
  p_depth

ggsave(
  filename = "merged_depth.pdf",
  plot = p_depth,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA/05-Liming/cellline/torun",
  width = 8,
  height = 5
)


# ADKP read_depth ---------------------------------------------------------

tibble::tibble(
  path = list.dirs(
    path = "/home/liuc9/github/scMOCHA/03-ADKP/output",
    full.names = T,
    recursive = F
  )
) |> 
  dplyr::filter(grepl(pattern = "R", x = path)) |> 
  dplyr::mutate(projectname = basename(path)) |> 
  dplyr::mutate(
    read_depth = purrr::map(
      .x = path,
      .f = \(.x) {
        # .x <- "/home/liuc9/github/scMOCHA/03-ADKP/output/R1246326"
        .read_depth <- data.table::fread(
          file = file.path(
            .x,
            "possorted_genome_bam.MT.depth"
          )
        )
      }
    )
  ) |> 
  dplyr::select(projectname, read_depth) |> 
  tidyr::unnest(cols = read_depth) ->
  adkp_read_depth


adkp_read_depth |> 
  dplyr::mutate(
    V3 = V3 / 1000000
  ) |> 
  ggplot(aes(
    x = V2,
    y = V3
  )) +
  geom_line(
    aes(color = projectname),
    stat = "identity"
  ) +
  geom_vline(xintercept = 3243, color = "black") +
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0)),
    limits = c(1, 17000),
    breaks = seq(0, 17000, 2000),
    labels = seq(0, 17000, 2000)
  ) +
  scale_y_continuous(
    expand = c(0.01, 0)
  ) +
  ggsci::scale_color_aaas(
    name = "Dataset",
    guide = guide_legend(nrow = 2)
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top"
  ) +
  labs(y = "Depth (x10^6)") ->
  adkp_p_depth;adkp_p_depth



tibble::tibble(
  path = list.dirs(
    path = "/home/liuc9/github/scMOCHA/01-Sci_Immunol_32651212/outputs",
    full.names = T,
    recursive = F
  )
) |> 
  # dplyr::filter(grepl(pattern = "R", x = path)) |> 
  dplyr::mutate(projectname = basename(path)) |> 
  dplyr::mutate(
    read_depth = purrr::map(
      .x = path,
      .f = \(.x) {
        # .x <- "/home/liuc9/github/scMOCHA/03-ADKP/output/R1246326"
        .read_depth <- data.table::fread(
          file = file.path(
            .x,
            "possorted_genome_bam.MT.depth"
          )
        )
      }
    )
  ) |> 
  dplyr::select(projectname, read_depth) |> 
  tidyr::unnest(cols = read_depth) ->
  sci_immunol_read_depth




sci_immunol_read_depth |> 
  dplyr::mutate(
    V3 = V3 / 1000000
  ) |> 
  ggplot(aes(
    x = V2,
    y = V3
  )) +
  geom_line(
    aes(color = projectname),
    stat = "identity"
  ) +
  geom_vline(xintercept = 3243, color = "black") +
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0)),
    limits = c(1, 17000),
    breaks = seq(0, 17000, 2000),
    labels = seq(0, 17000, 2000)
  ) +
  scale_y_continuous(
    expand = c(0.01, 0)
  ) +
  ggsci::scale_color_aaas(
    name = "Dataset",
    guide = guide_legend(nrow = 2)
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top"
  ) +
  labs(y = "Depth (x10^6)") ->
  sci_immunol_p_depth;sci_immunol_p_depth


# Mixed cellline ----------------------------------------------------------


gex_wt <- data.table::fread("/home/liuc9/github/scMOCHA/05-Liming/scmocha-celline/cromwell-executions/scMOCHA/1165f242-d68f-4260-8be4-55785cf2bc71/call-cellranger_count/execution/GEX_WT/outs/possorted_genome_bam.MT.depth") |> 
  tibble::as_tibble() |> 
  tibble::add_column(
    projectname = "GEX_WT",
    .before = 1
  )
# gex_rh0 <- data.table::fread("05-Liming/scmocha-celline/cromwell-executions/scMOCHA/4adb4525-bf49-445d-b433-40c4b84a9932/call-cellranger_count/execution/GEX_Rh0/outs/possorted_genome_bam.MT.depth") |> 
#   tibble::as_tibble() |> 
#   tibble::add_column(
#     projectname = "GEX_Rh0",
#     .before = 1
#   )
gex_read_depth <- gex_wt


gex_read_depth |> 
  dplyr::mutate(
    V3 = V3 / 1000000
  ) |> 
  ggplot(aes(
    x = V2,
    y = V3
  )) +
  geom_line(
    aes(color = projectname),
    stat = "identity"
  ) +
  geom_vline(xintercept = 3243, color = "black") +
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0)),
    limits = c(1, 17000),
    breaks = seq(0, 17000, 2000),
    labels = seq(0, 17000, 2000)
  ) +
  scale_y_continuous(
    expand = c(0.01, 0)
  ) +
  ggsci::scale_color_aaas(
    name = "Dataset",
    guide = guide_legend(nrow = 2)
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top"
  ) +
  labs(y = "Depth (x10^6)") ->
  gex_p_depth;gex_p_depth



# Combined all read depth -------------------------------------------------
dplyr::bind_rows(
  sci_immunol_read_depth |> 
    dplyr::group_by(V2) |> 
    dplyr::summarise(m = mean(V3)) |> 
    dplyr::mutate(
      project = "Sci_Immunol"
    ),
  gex_read_depth |> 
    dplyr::group_by(V2) |> 
    dplyr::summarise(m = mean(V3)) |> 
    dplyr::mutate(
      project = "Mixed 4 celllines WT"
    ),
  adkp_read_depth |> 
    dplyr::group_by(V2) |> 
    dplyr::summarise(m = mean(V3)) |> 
    dplyr::mutate(
      project = "ADKP"
    ),
  read_depth |> 
    dplyr::group_by(V2) |> 
    dplyr::summarise(m = mean(V3)) |> 
    dplyr::mutate(
      project = "Pei 143b"
    )
) ->
  merged_read_depth

merged_read_depth |> 
  dplyr::mutate(
    m = m / 1000000
  ) |>
  ggplot(aes(
    x = V2,
    y = m
  )) +
  geom_line(
    aes(color = project),
    stat = "identity"
  ) +
  geom_vline(xintercept = 3243, color = "black") +
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0)),
    limits = c(1, 17000),
    breaks = seq(0, 17000, 2000),
    labels = seq(0, 17000, 2000)
  ) +
  scale_y_continuous(
    expand = c(0.01, 0)
  ) +
  ggsci::scale_color_aaas(
    name = "Dataset"
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = c(0.8, 0.8),
    legend.key = element_blank()
  ) + 
  labs(y = "Depth (x10^6)") ->
  merged_p_read_depth;merged_p_read_depth


ggsave(
  filename = "combined_read_depth.pdf",
  plot = merged_p_read_depth,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA/05-Liming/cellline/torun",
  width = 8,
  height = 5
)



merged_read_depth |> 
  dplyr::filter(
    V2 > 3200
  ) |> 
  dplyr::filter(V2 < 3300) |> 
  # dplyr::mutate(
  #   m = m / 1000000
  # ) |>
  ggplot(aes(
    x = V2,
    y = m
  )) +
  geom_line(
    aes(color = project),
    stat = "identity"
  ) +
  geom_vline(xintercept = 3243, color = "black") +
  # scale_x_continuous(
  #   expand = expansion(mult = c(0.01, 0)),
  #   limits = c(1, 17000),
  #   breaks = seq(0, 17000, 2000),
  #   labels = seq(0, 17000, 2000)
  # ) +
  scale_y_continuous(
    expand = c(0.01, 0)
  ) +
  ggsci::scale_color_aaas(
    name = "Dataset"
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = c(0.8, 0.8),
    legend.key = element_blank()
  ) + 
  labs(y = "Depth (x10^6)") ->
  zoomin_merged_p_read_depth;zoomin_merged_p_read_depth


ggsave(
  filename = "zoomin_combined_read_depth.pdf",
  plot = zoomin_merged_p_read_depth,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA/05-Liming/cellline/torun",
  width = 8,
  height = 5
)


# Load mutations ----------------------------------------------------------

# 
# future::plan(future::multisession, workers = 7)
# outdir |> 
#   dplyr::mutate(
#     cluster = furrr::future_map(
#       .x = outputdir,
#       .f = \(.x) {
#         cell_variant_annotation <- readr::read_tsv(
#           file.path(
#             .x,
#             "cell_variant_annotation.tsv"
#           )
#         )
# 
#         cell_heteroplasmic_df_raw <- readr::read_tsv(
#           file.path(
#             .x,
#             "cell.cell_heteroplasmic_df_raw.tsv.gz"
#           )
#         )
#         cell_coverage <- readr::read_tsv(
#           file.path(
#             .x,
#             "cell.coverage.txt.gz"
#           )
#         )
#         
#         tibble::tibble(
#             cell_variant_annotation = list(cell_variant_annotation),
#             cell_heteroplasmic_df_raw = list(cell_heteroplasmic_df_raw),
#             cell_coverage = list(cell_coverage)
#         )
#       }
#     )
#   )  ->
#   mutations
# future::plan(future::sequential)
# 
# 
# mutations |> 
#   dplyr::select(projectname, cluster) |> 
#   tidyr::unnest(cols = cluster)

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(
  "/home/liuc9/github/scMOCHA/05-Liming/cellline/torun/05-stat-cellline.rda"
)
load("/home/liuc9/github/scMOCHA/05-Liming/cellline/torun/05-stat-cellline.rda")