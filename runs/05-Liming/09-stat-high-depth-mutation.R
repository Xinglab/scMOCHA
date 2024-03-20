#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Mar 19 17:38:14 2024
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
outdir <- readr::read_rds("/home/liuc9/github/scMOCHA/05-Liming/cellline/torun/outdir.rds.gz")

wt <- tibble::tibble(
  projectname = "WT",
  outputdir = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth/cromwell-executions/scMOCHA/5e46ec20-206d-443f-a390-e6507df10373/call-gather_outputfiles/execution/WT"
)
outdir <- dplyr::bind_rows(
  outdir, 
  wt
)

cellname <- c(
  "cluster_0" = "WAL2A",
  "cluster_1" = "A549",
  "cluster_2" = "HEK293",
  "cluster_3" = "143B"
)


# Variant intersection ----------------------------------------------------


future::plan(future::multisession, workers = 8)
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
        
        cell_heteroplasmic_df_raw <- data.table::fread(
          file.path(
            .x,
            "cluster.cell_heteroplasmic_df.tsv.gz"
          )
        )
        # cell_coverage <- data.table::fread(
        #   file.path(
        #     .x,
        #     "cluster.coverage.txt.gz"
        #   )
        # )
        
        tibble::tibble(
          cell_variant_annotation = list(cell_variant_annotation),
          cell_heteroplasmic_df_raw = list(cell_heteroplasmic_df_raw),
          # cell_coverage = list(cell_coverage)
        )
      }
    )
  )  ->
  mutations
future::plan(future::sequential)

mutations |> 
  dplyr::select(
    projectname,
    cluster
  ) |> 
  tidyr::unnest(cols = cluster) |> 
  dplyr::mutate(
    cell_variant_annotation_new = purrr::pmap(
      .l = list(
        projectname,
        cell_variant_annotation,
        cell_heteroplasmic_df_raw
      ),
      .f = \(projectname, cell_variant_annotation, cell_heteroplasmic_df_raw) {
        cell_variant_annotation |> 
          dplyr::mutate(v = "{Position}{Ref}>{Alt}" |> glue::glue()) ->
          cell_variant_annotation
        
        if(projectname != "WT") {return(
          cell_variant_annotation |> 
            dplyr::mutate(Sample = projectname) |> 
            dplyr::select(Sample,  v)
        )}
        
        cell_heteroplasmic_df_raw |>
          tidyr::gather(-V1, key = v, value = freq) |>
          dplyr::filter(freq > 0.05) ->
          .vv

        cell_variant_annotation |>
          dplyr::inner_join(.vv, by = "v") |>
          dplyr::select(Sample = V1, v)
      }
      
    )
  ) ->
  mutations_unnest

mutations_unnest |> 
  dplyr::select(cell_variant_annotation_new) |> 
  tidyr::unnest(cell_variant_annotation_new) |> 
  dplyr::rename(projectname = Sample) |> 
  dplyr::mutate(projectname = plyr::revalue(projectname, cellname)) |> 
  dplyr::group_by(v) |> 
  tidyr::nest() |> 
  dplyr::ungroup() |> 
  dplyr::mutate(
    srrid = purrr::map(
      .x = data,
      .f = \(.x) {
        .x |> dplyr::pull(projectname) ->
          .xx
      factor(
        .xx,
        levels = c(cellname, paste("Pei", 1:7, sep = "-"))
      )
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
  scale_x_upset(
    order_by = "degree",
    set = c(cellname, paste("Pei", 1:7, sep = "-"))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1), add = 0),
  ) +
  # axis_combmatrix()
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
ggplot2::ggsave(
  filename = "intersected_mutations.pdf",
  plot = intersected_mutations,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth",
  width = 12,
  height = 7
)



# single cell 3243 --------------------------------------------------------



future::plan(future::multisession, workers = 8)
outdir |>
  dplyr::mutate(
    n_mutations = furrr::future_map(
      .x = outputdir,
      .f = \(.x) {
        
        # .x <- outdir$outputdir[[8]]
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


barcode <- readr::read_tsv("/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth/cromwell-executions/scMOCHA/5e46ec20-206d-443f-a390-e6507df10373/call-gather_outputfiles/execution/WT/barcode_cluster.tsv", col_names = F) |> 
  dplyr::select(barcode = X1, X3) |> 
  dplyr::mutate(cellname = plyr::revalue(X3, cellname)) |> 
  dplyr::mutate(barcode = glue::glue("WT-{barcode}")) |> 
  dplyr::select(-X3)

agct_3243_ratio_cell |> 
  dplyr::mutate(barcode = glue::glue("{projectname}-{V2}")) |> 
  dplyr::filter(v %in% c("a", "g"))  |> 
  dplyr::select(-V2) |> 
  dplyr::left_join(barcode, by = "barcode") |> 
  dplyr::mutate(cellname = ifelse(is.na(cellname), projectname, cellname)) |> 
  dplyr::filter(cellname != "WT") |> 
  dplyr::rename(dataset = cellname) |> 
  dplyr::mutate(dataset = factor(dataset, levels =  c(paste("Pei", 1:7, sep = "-"), cellname |> rev()))) ->
  for_celllevel_plot

for_celllevel_plot |> 
  dplyr::select(-n) |> 
  tidyr::spread(key = v, value = ratio) |> 
  dplyr::arrange(dataset, -a, -g) ->
  sortv2

for_celllevel_plot |> 
  dplyr::mutate(
    barcode = factor(barcode, levels = sortv2$barcode)
  ) |> 
  ggplot(aes(
    x = barcode,
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
    barcode = factor(barcode, levels = sortv2$barcode)
  ) |> 
  dplyr::arrange(barcode) |> 
  dplyr::select(dataset, barcode) |> 
  dplyr::distinct() |> 
  tibble::rowid_to_column() |> 
  dplyr::group_by(dataset) |>
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
        string = dataset,
        width = 10
      ),
      ""
    )
  ) |> 
  ggplot(aes(x = barcode, y = 1, fill = dataset, label = label)) +
  geom_tile() +
  geom_text(
    color = "white",
    # angle = 90,
    size = 6
  ) +
  # ggsci::scale_fill_npg() 
  scale_fill_brewer(palette =  "Set3") +
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
  path = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth",
  width = 10,
  height = 5
)

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------