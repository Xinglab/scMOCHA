#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Mon Mar 18 14:47:46 2024
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
library(data.table)
#library(rlang)
library(GetoptLong)
library(logger)
library(gggenes)

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


# load data ---------------------------------------------------------------


mt_exons_df <- "/scr1/users/liuc9/mitochondrial/realdata/05-Liming/scmocha-mixed-cellline-high-depth/cromwell-executions/scMOCHA/023d7328-9097-4e50-8c11-19f860c5519e/call-cellranger_count/inputs/2014965526/mt_exons.df.rds.gz"
gtf_gene_df <-
  readr::read_rds(
    file = mt_exons_df
  )

depthfile <- "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth/cromwell-executions/scMOCHA/5e46ec20-206d-443f-a390-e6507df10373/call-gather_outputfiles/execution/WT/possorted_genome_bam.MT.depth"

coverage <- data.table::fread(
  input = depthfile,
  col.names = c("chr", "pos", "depth")
)

# body --------------------------------------------------------------------



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

coverage %>%
  dplyr::mutate(depth = depth / 10 ^6) |> 
  ggplot(aes(x = pos, y = depth)) +
  geom_bar(stat = "identity") +
  # geom_line() +
  # geom_vline(xintercept = 3243) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0)),
    limits = c(1, 17000),
  ) +
  scale_y_continuous(
    expand = c(0.01, 0), 
    limits = c(0, 2),
    # labels = scales::label_scientific(digits = 3)
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  ) +
  labs(y = latex2exp::TeX("$Depth (10^6$)") ) ->
  p1;p1

wrap_plots(
  p1,
  pg,
  ncol = 1,
  heights = c(0.9, 0.1)
) ->
  p;p

ggplot2::ggsave(
  filename = "coverage.pdf",
  plot = p,
  path = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth",
  width = 15,
  height = 7
)
gex_read_depth <- coverage |> 
  tibble::as_tibble() |> 
  tibble::add_column(
    projectname = "GEX_WT",
    .before = 1
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



# Sci_Immunol -------------------------------------------------------------


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

# Pei cell line -----------------------------------------------------------

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
gex_read_depth_old <- gex_wt

pei_read_depth <- readr::read_rds("/home/liuc9/github/scMOCHA/05-Liming/cellline/torun/pei-cellline-depth.rds.gz")

# Combined all read depth -------------------------------------------------
dplyr::bind_rows(
  sci_immunol_read_depth |> 
    dplyr::group_by(V2) |> 
    dplyr::summarise(m = mean(V3)) |> 
    dplyr::mutate(
      project = "Sci Immunol"
    ),
  gex_read_depth |> 
    dplyr::rename(V2 = pos, V3 = depth) |> 
    dplyr::group_by(V2) |> 
    dplyr::summarise(m = mean(V3)) |> 
    dplyr::mutate(
      project = "High depth mixed 4 celllines WT"
    ),
  gex_read_depth_old |> 
    dplyr::group_by(V2) |> 
    dplyr::summarise(m = mean(V3)) |> 
    dplyr::mutate(
      project = "Low depth mixed 4 celllines WT"
    ),
  adkp_read_depth |> 
    dplyr::group_by(V2) |> 
    dplyr::summarise(m = mean(V3)) |> 
    dplyr::mutate(
      project = "ADKP"
    ),
  pei_read_depth |> 
    dplyr::group_by(V2) |> 
    dplyr::summarise(m = mean(V3)) |> 
    dplyr::mutate(
      project = "Pei 143B"
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
    expand = c(0.01, 0), 
    limits = c(0, 2),
    # labels = scales::label_scientific(digits = 3)
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
    legend.position = c(0.8, 0.8),
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
  labs(y = latex2exp::TeX("$Depth (10^6$)") ) ->
  merged_p_read_depth;merged_p_read_depth

wrap_plots(
  merged_p_read_depth,
  pg,
  ncol = 1,
  heights = c(0.9, 0.1)
) ->
  pg_merged_p_read_depth;pg_merged_p_read_depth

ggplot2::ggsave(
  filename = "combined_read_depth.pdf",
  plot = pg_merged_p_read_depth,
  path = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth",
  width = 15,
  height = 7
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
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.line.x = element_blank(),
    # axis.title.x = element_blank(),
    legend.position = c(0.8, 0.8),
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
    ),
    # axis.line.y.left = element_line(color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.line.x = element_blank(),
    axis.title.x = element_blank(),
  ) + 
  labs(y = latex2exp::TeX("$Depth (10^6$)") )  ->
  zoomin_merged_p_read_depth;zoomin_merged_p_read_depth
ggplot2::ggsave(
  filename = "zoomin_combined_read_depth.pdf",
  plot = zoomin_merged_p_read_depth,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth",
  width = 12,
  height = 5
)

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image("/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth/08-stat-high-depth-scRNA-read-coverage.rda")
