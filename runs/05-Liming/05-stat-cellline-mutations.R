#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Dec 14 13:51:42 2023
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

outdir <- readr::read_rds("/home/liuc9/github/scMOCHA/05-Liming/cellline/torun/outdir.rds.gz")
# body --------------------------------------------------------------------



future::plan(future::multisession, workers = 7)
outdir |>
  dplyr::mutate(
    cluster = furrr::future_map(
      .x = outputdir,
      .f = \(.x) {
        cell_variant_annotation <- data.table::fread(
          file.path(
            .x,
            "cell_variant_annotation.tsv"
          )
        )

        # cell_heteroplasmic_df_raw <- data.table::fread(
        #   file.path(
        #     .x,
        #     "cell.cell_heteroplasmic_df_raw.tsv.gz"
        #   )
        # )
        # cell_coverage <- data.table::fread(
        #   file.path(
        #     .x,
        #     "cell.coverage.txt.gz"
        #   )
        # )

        tibble::tibble(
            cell_variant_annotation = list(cell_variant_annotation),
            # cell_heteroplasmic_df_raw = list(cell_heteroplasmic_df_raw),
            # cell_coverage = list(cell_coverage)
        )
      }
    )
  )  ->
  mutations
future::plan(future::sequential)


# Up set plot -------------------------------------------------------------


mutations |> 
  dplyr::select(
    projectname,
    cluster
  ) |> 
  tidyr::unnest(cols = cluster) ->
  mutations_unnest

mutations_unnest |> 
  dplyr::select(projectname, cell_variant_annotation) |> 
  tidyr::unnest(cols = cell_variant_annotation) |> 
  dplyr::mutate(v = "{Position}{Ref}>{Alt}" |> glue::glue()) ->
  mutations_unnest_v

mutations_unnest_v |> 
  dplyr::select(projectname, v) |> 
  dplyr::group_by(v) |> 
  tidyr::nest() |> 
  dplyr::ungroup() |> 
  dplyr::mutate(
    srrid = purrr::map(
      .x = data,
      .f = \(.x) {
        .x |> dplyr::pull(projectname)
      }
    )
  ) |> 
  dplyr::select(-data) ->
  for_upset

library(ggupset)
for_upset |> 
  ggplot(aes(x = srrid)) +
  geom_bar(width = 0.5, fill = ggsci::pal_jama()(4)[[1]]) +
  geom_text(
    stat='count',
    aes(label=after_stat(count)),
    vjust = -0.5,
    color = "black",
    size = 6,
    fontface = "bold"
  ) +
  scale_x_upset(order_by = "degree") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1), add = 0)
  ) +
  theme_combmatrix(
    combmatrix.label.make_space = TRUE,
    combmatrix.panel.point.color.fill =  ggsci::pal_jama()(4)[[1]],
    combmatrix.panel.line.size = 0,
    combmatrix.label.text = element_text(
      size = 12,
      color = "black",
      face = "bold"
    ),
    combmatrix.label.extra_spacing = 5,
    combmatrix.panel.striped_background.color.one = "white",
    combmatrix.panel.striped_background.color.two = "grey",
  ) +
  labs(
    y = "# of Variants",
    x = "",
    # title = .x
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.title.y = element_text(
      size = 16,
      color = "black",
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 14,
      color = "black"
    ),
    plot.title = element_text(
      hjust = 0.5,
      color = "black",
      size = 16,
      face = "bold"
    )
  ) ->
  intersected_mutations;intersected_mutations


ggsave(
  filename = "intersected_mutations.pdf",
  plot = intersected_mutations,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA/05-Liming/cellline/torun",
  width = 7,
  height = 5
)


# 3243 --------------------------------------------------------------------



future::plan(future::multisession, workers = 7)
outdir |>
  dplyr::mutate(
    n_mutations = furrr::future_map(
      .x = outputdir,
      .f = \(.x) {
        A <- data.table::fread(
          file.path(
            .x,
            "cell.A.txt.gz"
          )
        ) |> 
          dplyr::filter(V1 == 3243)
        n_A <- sum(A$V3, A$V4)
        
        G <- data.table::fread(
          file.path(
            .x,
            "cell.G.txt.gz"
          )
        ) |> 
          dplyr::filter(V1 == 3243)
        n_G <- sum(G$V3, G$V4)
        
        C <- data.table::fread(
          file.path(
            .x,
            "cell.C.txt.gz"
          )
        ) |> 
          dplyr::filter(V1 == 3243)
        
        n_C <- sum(C$V3, C$V4)
        
        TT <- data.table::fread(
          file.path(
            .x,
            "cell.T.txt.gz"
          )
        ) |> 
          dplyr::filter(V1 == 3243)
        n_T <- sum(TT$V3, TT$V4)
        
        tibble::tibble(
          a = n_A,
          g = n_G,
          c = n_C,
          t = n_T
        )
      }
    )
  )  ->
  agct_3243
future::plan(future::sequential)

agct_3243 |> 
  dplyr::select(projectname, n_mutations) |> 
  tidyr::unnest(n_mutations) |> 
  tidyr::pivot_longer(
    cols = -projectname,
    names_to = "v",
    values_to = "n"
  ) |> 
  dplyr::filter(
    v %in% c("a", "g")
  ) |> 
  dplyr::group_by(
    projectname
  ) |> 
  dplyr::mutate(
    ratio = n / sum(n) 
  ) |> 
  dplyr::ungroup() ->
  agct_3243_ratio

agct_3243_ratio |> 
  dplyr::filter(projectname == "Pei-5")

agct_3243_ratio |> 
  dplyr::filter(v %in% c("a", "g")) |> 
  dplyr::select(-ratio) |> 
  tidyr::pivot_wider(
    names_from = v,
    values_from = n
  ) |> 
  dplyr::mutate(label = "A({round(a / (a + g) * 100,2)}%)\nG({round(g / (a + g) * 100, 2)}%)\n(n={a+g})" |> glue::glue()) ->
  forlabel

agct_3243_ratio |> 
  dplyr::filter(v %in% c("a", "g")) |> 
  ggplot(aes(
    x = projectname,
    y = ratio
  )) +
  geom_col(aes(fill = v)) +
  ggsci::scale_fill_aaas(
    name = "Variant",
    label = c("A", "G")
  ) +
  geom_text(
    data = forlabel,
    aes(
      x = projectname,
      y = 1,
      label = label
    ),
    color = "black",
    vjust = 0,
    size = 4
  ) +
  scale_y_continuous(
    expand = expansion(mult = 0, add = c(0.01, 0.2)),
    breaks = seq(0, 1, 0.2),
    labels = seq(0, 1, 0.2)
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    # axis.ticks.x = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 14
    ),
    axis.title.y = element_text(
      color = "black",
      size = 16
    ),
    # axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.key = element_blank()
  ) +
  labs(
    y = "Ratio"
  ) ->
  ag_ratio;ag_ratio

ggsave(
  filename = "ag_ratio.pdf",
  plot = ag_ratio,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA/05-Liming/cellline/torun",
  width = 7,
  height = 5
)


# Single cell level -------------------------------------------------------


future::plan(future::multisession, workers = 7)
outdir |>
  dplyr::mutate(
    n_mutations = furrr::future_map(
      .x = outputdir,
      .f = \(.x) {
        
        # .x <- outdir$outputdir[[7]]
        A <- data.table::fread(
          file.path(
            .x,
            "cell.A.txt.gz"
          )
        ) |> 
          dplyr::filter(V1 == 3243) |> 
          dplyr::mutate(a = V3 + V4) |> 
          dplyr::select(V2, a)
        
        G <- data.table::fread(
          file.path(
            .x,
            "cell.G.txt.gz"
          )
        ) |> 
          dplyr::filter(V1 == 3243) |> 
          dplyr::mutate(g = V3 + V4) |> 
          dplyr::select(V2, g)
        
        C <- data.table::fread(
          file.path(
            .x,
            "cell.C.txt.gz"
          )
        ) |> 
          dplyr::filter(V1 == 3243) |> 
          dplyr::mutate(c = V3 + V4) |> 
          dplyr::select(V2, c)
        
        
        TT <- data.table::fread(
          file.path(
            .x,
            "cell.T.txt.gz"
          )
        ) |> 
          dplyr::filter(V1 == 3243) |> 
          dplyr::mutate(t = V3 + V4) |> 
          dplyr::select(V2, t)
        
        dplyr::full_join(
          A, G, by = "V2"
        ) |> 
          dplyr::full_join(
            C, by = "V2"
          ) |> 
          dplyr::full_join(
            TT, by = "V2"
          ) |> 
          tidyr::pivot_longer(
            cols = -V2,
            names_to = "v",
            values_to = "n"
          ) |> 
          tidyr::replace_na(
            replace = list(n = 0)
          ) |> 
          dplyr::group_by(V2) |> 
          dplyr::mutate(
            ratio = n / sum(n)
          ) |> 
          dplyr::ungroup()
      }
    )
  )  ->
  agct_3243_cell
future::plan(future::sequential)

agct_3243_cell |> 
  dplyr::select(projectname, n_mutations) |> 
  tidyr::unnest(n_mutations) ->
  agct_3243_ratio_cell


agct_3243_ratio_cell |> 
  dplyr::mutate(V2 = glue::glue("{projectname}-{V2}")) |> 
  dplyr::filter(v %in% c("a", "g"))  ->
  for_celllevel_plot

for_celllevel_plot |> 
  dplyr::select(-n) |> 
  tidyr::spread(key = v, value = ratio) |> 
  dplyr::arrange(projectname, -a, -g) ->
  sortv2

for_celllevel_plot |> 
  dplyr::mutate(
    V2 = factor(V2, levels = sortv2$V2)
  ) |> 
  ggplot(aes(
    x = V2,
    y = ratio
  )) +
  geom_col(aes(fill = v))  +
  ggsci::scale_fill_aaas(
    name = "Variant",
    label = c("A", "G")
  ) +
  scale_y_continuous(
    expand = expansion(mult = 0, add = c(0.01, 0.01)),
    breaks = seq(0, 1, 0.2),
    labels = seq(0, 1, 0.2)
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    # axis.line.x.bottom = element_line(color = "black"),
    # axis.ticks.x = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 14
    ),
    axis.title.y = element_text(
      color = "black",
      size = 16
    ),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.key = element_blank(),
  ) +
  labs(
    y = "Ratio"
  ) ->
  p_col;p_col

for_celllevel_plot |> 
  dplyr::mutate(
    V2 = factor(V2, levels = sortv2$V2)
  ) |> 
  dplyr::arrange(V2) |> 
  dplyr::select(projectname, V2) |> 
  dplyr::distinct() |> 
  tibble::rowid_to_column() |> 
  dplyr::group_by(projectname) |>
  dplyr::mutate(
    n = mean(rowid)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    n = as.integer(floor(n))
  ) |> 
  dplyr::mutate(
    label = ifelse(
      rowid == n,
      stringr::str_wrap(
        string = projectname,
        width = 10
      ),
      ""
    )
  ) |> 
  ggplot(aes(x = V2, y = 1, fill = projectname, label = label)) +
  geom_tile() +
  geom_text(
    color = "white",
    # angle = 90,
    size = 6
  ) +
  ggsci::scale_fill_npg() +
  scale_y_continuous(
    expand = c(0,0),
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
  ) ->
  p_bar;p_bar
  


p_col/  p_bar + plot_layout(heights = c(9, 1)) ->
  p_bar_col;p_bar_col

ggsave(
  filename = "cell_ag_ratio.pdf",
  plot = p_bar_col,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA/05-Liming/cellline/torun",
  width = 10,
  height = 5
)


# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(
  "/home/liuc9/github/scMOCHA/05-Liming/cellline/torun/06-stat-cellline-mutations.rda"
)
