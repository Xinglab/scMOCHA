#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Fri Aug 18 18:11:28 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
library(Seurat)

# library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------

outfiles <- readr::read_tsv(
  file = "/home/liuc9/github/scMOCHA/03-ADKP/output/outfiles.tsv"
)

outfiles$outfile[[2]]

supdata1_1 <- readxl::read_xls(
  path = "/scr1/users/liuc9/mitochondrial/realdata/03-ADKP/output/OlhaEtAl-NatCommu-2020-supdata1.xls",
  sheet = 1,
  n_max = 18
) |>
  dplyr::mutate(
    `Sample ID` = gsub(
      pattern = " GM",
      replacement = "",
      x = `Sample ID`
    )
  )


supdata1_2 <- readxl::read_xls(
  path = "/scr1/users/liuc9/mitochondrial/realdata/03-ADKP/output/OlhaEtAl-NatCommu-2020-supdata1.xls",
  sheet = 2,
  n_max = 15
)

synapse_metadata <- readr::read_tsv(
  file = "/home/liuc9/github/scMOCHA/03-ADKP/SYNAPSE_METADATA_MANIFEST.tsv"
)

synapse_metadata |>
  dplyr::glimpse()

synapse_metadata |>
  dplyr::select(
    specimenID,
    sex,
    cellType,
    diagnosis,
    individualID,
  ) |>
  dplyr::distinct() |>
  dplyr::arrange(specimenID, individualID) |>
  dplyr::mutate(
    `Sample ID` = gsub(
      pattern = "Microglia_MO_",
      replacement = "",
      x = specimenID
    )
  ) |>
  dplyr::mutate(
    srrid = individualID
  ) ->
srarun

outfiles |>
  dplyr::left_join(
    srarun,
    by = "srrid"
  ) |>
  dplyr::left_join(
    supdata1_1,
    by = "Sample ID"
  ) |>
  dplyr::left_join(
    supdata1_2,
    by = "Sample ID"
  ) ->
outfiles_supdata




# body --------------------------------------------------------------------


outfiles_supdata |>
  dplyr::mutate(
    a = purrr::map(
      .x = tardir,
      .f = \(.x) {
        if (is.na(.x)) {
          return(NA)
        }

        readxl::read_excel(
          path = file.path(
            .x, "qc_cell_stats.xlsx"
          )
        )
      }
    )
  ) |>
  tidyr::unnest(cols = a) ->
metadata


readr::write_csv(
  x = metadata,
  file = "/home/liuc9/github/scMOCHA/03-ADKP/output/metadata.csv"
)




# Load mutation --------------------------------------------------------------

metadata |>
  dplyr::mutate(
    anno = purrr::map(
      .x = outfile,
      .f = \(.x) {
        if (.x == "FALSE") {
          return(NULL)
        }

        .uuid <- dirname(dirname(dirname(.x)))

        .cva <- file.path(
          .uuid,
          "call-plot_scMOCHA/execution",
          "cell_variant_annotation.tsv"
        )

        readr::read_tsv(
          file = .cva,
          show_col_types = FALSE
        )
      }
    )
  ) |>
  dplyr::mutate(
    nmut = purrr::map_int(
      .x = anno,
      .f = \(.x) {
        if (is.null(.x)) {
          return(NA_integer_)
        }
        nrow(.x)
      }
    )
  ) |>
  dplyr::mutate(
    haplogroup = purrr::map2(
      .x = anno,
      .y = srrid,
      .f = \(.x, .y) {
        message(.y)
        if (is.null(.x)) {
          return(
            tibble::tibble(
              Haplogroup = NA_character_,
              Verbose_haplogroup = NA_character_
            )
          )
        }
        .x |>
          dplyr::select(Haplogroup, Verbose_haplogroup) |>
          dplyr::filter(!is.na(Haplogroup)) |>
          dplyr::distinct() |>
          dplyr::mutate_all(.funs = as.character)
      }
    )
  ) |>
  tidyr::unnest(cols = haplogroup) ->
metadata_anno

readr::write_rds(
  x = metadata_anno,
  file = "/home/liuc9/github/scMOCHA/03-ADKP/output/metadata_anno.rds"
)

metadata_anno |> dplyr::glimpse()

# Select output columns ---------------------------------------------------


metadata_anno |>
  dplyr::mutate(
    pass = ifelse(
      test = is.na(tardir),
      yes = "Fail",
      no = "Pass"
    )
  ) |>
  dplyr::mutate(
    ratio = round(`number of cells after filtering` / `estimated number of cells`, 2)
  ) |>
  dplyr::select(
    srrid,
    `Sample ID`,
    Age, Sex,
    `Diagnosis (neurology)`,
    `Gating strategy`,
    study,
    `Median UMI/cell` = `median UMI counts per cell`,
    `Median genes/cell` = `median genes per cell`,
    `# of cells` = `estimated number of cells`,
    `# cells after filter` = `number of cells after filtering`,
    `Cell ratio` = ratio,
    `# of variants` = nmut,
    Haplogroup = Haplogroup,
    Haplogroup_v = Verbose_haplogroup
  ) |>
  dplyr::arrange(`Sample ID`) ->
metadata_clean

metadata_clean |>
  writexl::write_xlsx(
    path = "/home/liuc9/github/scMOCHA/03-ADKP/output/metadata_clean.xlsx"
  )


# Single cell consistency -------------------------------------------------

future::plan(future::multisession, workers = 10)
metadata_anno |>
  dplyr::mutate(
    sc_azimuth = furrr::future_map(
      .x = tardir,
      .f = \(.x) {
        if (is.na(.x)) {
          return(NULL)
        }

        .sc_a <- file.path(.x, "sc_azimuth.rds.gz") |>
          readr::read_rds()
        .sc_a
      }
    )
  ) ->
metadata_anno_azimuth
future::plan(future::sequential)

raw_annotations <- readr::read_csv(
  file = "/scr1/users/liuc9/mitochondrial/realdata/03-ADKP/forrefs/annofile.csv"
) |>
  dplyr::mutate(cluster = glue::glue("cluster_{cluster_label}")) |>
  dplyr::mutate(cluster = forcats::fct_reorder(cluster, cluster_label)) |>
  dplyr::select(cellid = sample_id, cluster)

cell_barcodes <- readr::read_csv(
  file = "/scr1/users/liuc9/mitochondrial/realdata/03-ADKP/forrefs/cell_barcodes.csv"
) |>
  dplyr::select(-1, -batchname) |>
  dplyr::left_join(raw_annotations, by = "cellid") |>
  dplyr::mutate(
    donor = gsub(
      pattern = "_GM",
      replacement = "",
      x = donor
    )
  ) |>
  dplyr::rename(`Sample ID` = donor) |>
  dplyr::group_by(`Sample ID`) |>
  tidyr::nest() |>
  dplyr::ungroup()

metadata_anno_azimuth |>
  dplyr::select(srrid, `Sample ID`, sc_azimuth) |>
  dplyr::left_join(cell_barcodes, by = "Sample ID") |>
  dplyr::mutate(
    a = purrr::map2(
      .x = sc_azimuth,
      .y = data,
      .f = \(.x, .y) {
        # .x <- d$sc_azimuth[[2]]
        # .y <- d$data[[2]]

        if (is.null(.x)) {
          return(
            tibble::tibble(
              n_raw = NULL,
              n_azi = NULL,
              n_nc = nrow(.y),
              azi_nc = NULL
            )
          )
        }

        .xx <- .x$cellbarcode_cluster |>
          dplyr::select(bc = cellbarcode, azi_cluster = cluster)

        n_raw <- nrow(.x$sc@meta.data)
        n_azi <- nrow(.xx)
        n_nc <- nrow(.y)

        .xx |>
          dplyr::left_join(.y, by = "bc") |>
          dplyr::mutate(
            azi_cluster = as.character(azi_cluster),
            cluster = as.character(cluster)
          ) ->
        .xxx

        tibble::tibble(
          n_raw = n_raw,
          n_azi = n_azi,
          n_nc = n_nc,
          azi_nc = list(.xxx)
        )
      }
    )
  ) |>
  tidyr::unnest(cols = a) ->
metadata_anno_azimuth_nc

metadata_anno_azimuth_nc |>
  dplyr::select(-sc_azimuth, -data) |>
  dplyr::mutate(
    ratio1 = round(n_azi / n_raw, 2),
    ratio2 = round(n_nc / n_raw, 2)
  ) |>
  dplyr::select(
    ID = srrid,
    `Sample ID`,
    `# Raw` = n_raw,
    `# scMOCHA` = n_azi,
    `# NC` = n_nc,
    `Ratio of scMOCHA` = ratio1,
    `Ratio of NC` = ratio2
  ) |>
  dplyr::arrange(`Sample ID`) |>
  writexl::write_xlsx(
    path = "/home/liuc9/github/scMOCHA/03-ADKP/output/cell_n_consistence.xlsx"
  )


# Cell consistence --------------------------------------------------------



metadata_anno_azimuth_nc |>
  dplyr::select(srrid, `Sample ID`, azi_nc) |>
  tidyr::unnest(cols = azi_nc) |>
  dplyr::filter(!is.na(cellid)) ->
metadata_anno_azimuth_nc_unnest

metadata_anno_azimuth_nc_unnest |>
  dplyr::count(azi_cluster, cluster) ->
forplot

forplot$n |> sum()

forplot |>
  ggplot(aes(
    x = azi_cluster,
    y = cluster,
  )) +
  geom_tile(aes(fill = n)) +
  scale_fill_gradient2(
    low = "grey",
    mid = "gold",
    high = "#F02415",
    midpoint = 500,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  ) +
  geom_text(
    aes(label = n),
    color = "white"
  ) +
  scale_x_discrete(
    limit = paste("cluster", c(1:11, 13), sep = "_"),
    expand = expansion(mult = 0, add = 0)
  ) +
  scale_y_discrete(
    limit = paste("cluster", c(1:11, 13), sep = "_"),
    expand = expansion(mult = 0, add = 0)
  ) +
  geom_segment(
    x = 1,
    y = 0,
    xend = 13,
    yend = 12,
    color = "grey"
  ) +
  geom_segment(
    x = 0,
    y = 1,
    xend = 12,
    yend = 13,
    color = "grey"
  ) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    # legend.position = "top",
    aspect.ratio = 1,
    axis.title = element_text(
      color = "black",
      size = 14
    ),
    plot.title = element_text(
      color = "black",
      size = 16
    )
  ) +
  labs(
    x = "scMOCHA cell cluster",
    y = "Olha et al cell cluster",
    title = "Total number of cells n = {sum(forplot$n)}" |> glue::glue()
  ) ->
tileplot

ggsave(
  filename = "cell-consistence-tileplot.pdf",
  plot = tileplot,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA/03-ADKP/output",
  width = 8,
  height = 7
)

# multi-class metrics -----------------------------------------------------


true_labels <- as.factor(metadata_anno_azimuth_nc_unnest$cluster)

predicted_labels <- factor(metadata_anno_azimuth_nc_unnest$azi_cluster, level = levels(true_labels))

classifier_metrics <- mltest::ml_test(predicted_labels, true_labels, output.as.table = FALSE)

# overall classification accuracy
accuracy <- classifier_metrics$accuracy

# F1-measures for classes "cat", "dog" and "rat"
F1 <- classifier_metrics$F1

# tabular view of the metrics (except for 'accuracy' and 'error.rate')
classifier_metrics <- mltest::ml_test(predicted_labels, true_labels, output.as.table = TRUE)



# Cell ratio --------------------------------------------------------------

metadata_anno |>
  dplyr::mutate(
    cellratio = purrr::map(
      .x = tardir,
      .f = function(.x) {
        if (is.na(.x)) {
          return(NULL)
        }
        .ratio <- readr::read_tsv(
          file = file.path(
            .x, "celltype_ratio.tsv"
          ),
          show_col_types = FALSE
        )
        .ratio
      }
    )
  ) ->
metadata_anno_cellratio

metadata_anno_cellratio |>
  dplyr::filter(!purrr::map_lgl(.x = cellratio, .f = is.null)) |>
  dplyr::select(`Sample ID`, dia = `Diagnosis (neurology)`, cellratio) |>
  dplyr::arrange(dplyr::desc(`Sample ID`)) ->
for_ratio_plot

for_ratio_plot |>
  tidyr::unnest(cellratio) |>
  dplyr::mutate(celltype = factor(
    celltype,
    levels = paste("cluster", c(1:11, 13), sep = "_")
  )) |>
  ggplot(aes(
    x = `ratio`,
    y = `Sample ID`,
    fill = celltype
  )) +
  geom_col() +
  scale_fill_manual(
    values = paletteer::paletteer_d(
      palette = "ggsci::springfield_simpsons",
      direction = -1
    ),
    # limits = paste("cluster", c(1:11, 13), sep = "_"),
    name = "Cell type"
  ) +
  scale_x_continuous(
    expand = expansion(mult = 0, add = 0),
    labels = scales::percent_format()
  ) +
  scale_y_discrete(
    limits = for_ratio_plot$`Sample ID`
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size = 12, face = "bold"),
    axis.text.y = element_text(
      color = for_ratio_plot$color
    ),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Cell ratio") ->
p_cellratio
p_cellratio

ggsave(
  filename = "Cell_ratio.pdf",
  plo = p_cellratio,
  device = "pdf",
  width = 11,
  height = 5,
  path = "/home/liuc9/github/scMOCHA/03-ADKP/output"
)


# Depth ----------------------------------------------------------------

gtf_gene_df <-
  readr::read_rds(
    file = "/home/liuc9/github/scMOCHA/fasta/mt_exons.df.rds.gz"
  )

# Check if ggtranscript is installed, install if not
if (!requireNamespace("ggtranscript", quietly = TRUE)) {
  message("Installing ggtranscript from GitHub...")
  devtools::install_github("dzhang32/ggtranscript")
}
library(ggtranscript)
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
p_mt_chrom
p_mt_chrom

metadata_anno_cellratio |>
  dplyr::mutate(
    depth = purrr::map(
      .x = tardir,
      .f = function(.x) {
        if (is.na(.x)) {
          return(NULL)
        }
        data.table::fread(
          input = file.path(
            .x, "possorted_genome_bam.MT.depth"
          ),
          col.names = c("chr", "pos", "depth")
        )
      }
    )
  ) ->
metadata_anno_depth

metadata_anno_depth |>
  dplyr::filter(!purrr::map_lgl(.x = depth, .f = is.null)) |>
  # dplyr::select(srrid, source_name, depth) |>
  dplyr::select(`Sample ID`, dia = `Diagnosis (neurology)`, depth) |>
  dplyr::arrange(dplyr::desc(`Sample ID`)) |>
  dplyr::mutate(`Sample ID` = factor(`Sample ID`)) |>
  dplyr::mutate(color = dplyr::case_match(
    dia,
    "MCI" ~ ggsci::pal_jama()(4)[[1]],
    "AD" ~ ggsci::pal_jama()(4)[[2]],
  )) ->
for_depth_plot


for_depth_plot |>
  tidyr::unnest(cols = depth) |>
  ggplot(aes(x = pos, y = depth, fill = `color`)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.03)),
    limits = c(1, 17000),
  ) +
  scale_y_continuous(
    expand = c(0.01, 0),
  ) +
  scale_fill_identity(
    name = "Sample"
  ) +
  # scale_fill_manual(
  #   name = "Sample",
  #   values = for_depth_plot$color,
  #   guide = guide_legend(nrow = 3)
  # )
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
    # strip.text = element_text(
    #   color = for_depth_plot$color
    # ),
    strip.text = element_blank(),
    legend.position = "none"
  ) +
  facet_wrap(
    facets = ~`Sample ID`,
    ncol = 1,
    strip.position = "right"
  ) +
  labs(y = "Depth") ->
p_mt_depth
p_mt_depth

p_depth <- cowplot::plot_grid(
  plotlist = list(p_mt_depth, p_mt_chrom),
  ncol = 1,
  align = "v",
  rel_heights = c(0.7, 0.3)
)

ggsave(
  filename = "Sample_depth_merge.pdf",
  plo = p_depth,
  device = "pdf",
  width = 15,
  height = 15,
  path = "/home/liuc9/github/scMOCHA/03-ADKP/output"
)

metadata_anno_depth$depth[[4]]$depth |> summary()

readr::write_rds(
  x = metadata_anno_depth,
  file = "/home/liuc9/github/scMOCHA/03-ADKP/output/rda/metadata_anno_depth.rds"
)

# Correlation 3D -------------------------------------------------------------



metadata_anno_depth$`estimated number of cells`
metadata_anno_depth$`number of cells after filtering`
metadata_anno_depth$`median UMI counts per cell`
metadata_anno_depth$nmut




metadata_anno_depth |>
  dplyr::filter(!is.na(nmut)) |>
  ggplot(aes(
    x = `number of cells after filtering`,
    y = nmut
  )) +
  geom_point()

metadata_anno_depth |>
  dplyr::filter(!is.na(nmut)) |>
  ggplot(aes(
    x = `median UMI counts per cell`,
    y = nmut
  )) +
  geom_point()

library(plotly)
plot_ly(
  data = metadata_anno_depth |>
    dplyr::mutate(dia = `Diagnosis (neurology)`) |>
    dplyr::filter(!is.na(nmut)),
  x = ~`median UMI counts per cell`,
  y = ~`number of cells after filtering`,
  z = ~nmut,
  # size = ~nmut,
  color = ~dia,
  colors = ggsci::pal_jama()(4)[c(2, 1)]
) |>
  layout(
    scene = list(
      xaxis = list(
        title = "Median UMI counts per cell",
        gridcolor = "rgb(255, 255, 255)",
        # range = c(2.003297660701705, 5.191505530708712),
        # type = 'log',
        zerolinewidth = 1,
        ticklen = 5,
        gridwidth = 2
      ),
      yaxis = list(
        title = "Number of cells after filtering",
        gridcolor = "rgb(255, 255, 255)",
        # range = c(36.12621671352166, 91.72921793264332),
        zerolinewidth = 1,
        ticklen = 5,
        gridwith = 2
      ),
      zaxis = list(
        title = "Number of Mutations",
        gridcolor = "rgb(255, 255, 255)",
        # type = 'log',
        zerolinewidth = 1,
        ticklen = 5,
        gridwith = 2
      )
    ),
    paper_bgcolor = "rgb(243, 243, 243)",
    plot_bgcolor = "rgb(243, 243, 243)"
  ) |>
  plotly::config(
    displayModeBar = TRUE,
    showEditInChartStudio = TRUE,
    plotlyServerURL = "https://chart-studio.plotly.com",
    displaylogo = FALSE
  ) ->
p3d
p3d

# reticulate::py_run_string("import sys")

htmlwidgets::saveWidget(
  p3d,
  file = file.path(
    "/home/liuc9/github/scMOCHA/03-ADKP/output",
    "nmut-vs-umi-ncell.html"
  ),
)



metadata_anno_depth$`number of cells after filtering`
metadata_anno_depth$`median UMI counts per cell`
metadata_anno_depth$nmut





# all factor correlations -------------------------------------------------


metadata_anno_depth |> colnames()
readr::write_rds(
  x = metadata_anno_depth,
  file = "/home/liuc9/github/scMOCHA/03-ADKP/output/rda/metadata_anno_depth.rds"
)

metadata_anno_depth |>
  dplyr::mutate(dia = `Diagnosis (neurology)`) |>
  dplyr::select(
    srrid,
    dia,
    Age,
    Sex,
    nmut,
    `median UMI counts per cell`,
    `number of cells after filtering`,
    depth
  ) |>
  dplyr::filter(!is.na(nmut)) |>
  dplyr::mutate(
    n_na = purrr::map(
      .x = depth,
      .f = \(.d) {
        .d |>
          dplyr::summarise(
            dep_s = sum(depth),
            dep_mea = mean(depth),
            dep_med = median(depth)
          )
      }
    )
  ) |>
  dplyr::select(-depth) |>
  tidyr::unnest(cols = n_na) |>
  dplyr::mutate(
    Sex = factor(Sex),
    dia = factor(dia)
  ) ->
metadata_anno_depth_dep



metadata_anno_depth_dep

t.test(nmut ~ dia, data = metadata_anno_depth_dep) |> report::report()
t.test(nmut ~ Sex, data = metadata_anno_depth_dep) |> report::report()

glm(
  nmut ~ dia + dep_med,
  data = metadata_anno_depth_dep
) |>
  report::report()

correlation::correlation(
  metadata_anno_depth_dep |>
    dplyr::select(-dep_s, -dep_mea),
  p_adjust = "none"
) |>
  summary(redundant = TRUE) ->
cor_summr

plot(cor_summr)

cor_summr |>
  as.data.frame() |>
  tidyr::pivot_longer(
    cols = -Parameter,
    names_to = "var2",
    values_to = "pval"
  ) |>
  dplyr::filter(!is.na(pval)) |>
  dplyr::filter(Parameter != var2) |>
  ggplot(aes(
    x = Parameter,
    y = var2,
    fill = pval
  )) +
  geom_tile() +
  geom_text(aes(label = round(pval, 2))) +
  scale_fill_gradient2(
    name = "R",
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    # space = "Lab",
    # na.value = "grey50",
    guide = "colourbar",
    # aesthetics = "colour"
  ) +
  scale_x_discrete(
    limits = c("nmut", "dep_med", "number of cells after filtering", "median UMI counts per cell", "Age"),
    labels = c("# variants", "median depth", "# cells", "median UMI/cell", "Age") |> stringr::str_to_sentence()
  ) +
  scale_y_discrete(
    limits = c("nmut", "dep_med", "number of cells after filtering", "median UMI counts per cell", "Age") |> rev(),
    labels = c("# variants", "median depth", "# cells", "median UMI/cell", "Age") |> rev() |> stringr::str_to_sentence()
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 14
    )
  ) ->
p_cor
ggsave(
  filename = "All-factor-correlations.pdf",
  plo = p_cor,
  device = "pdf",
  width = 9,
  height = 6,
  path = "/home/liuc9/github/scMOCHA/03-ADKP/output"
)


cor.test(
  formula = ~ nmut + dep_med,
  data = metadata_anno_depth_dep
) ->
ct


metadata_anno_depth_dep |>
  # dplyr::filter(dep_med > 1000) |>
  dplyr::mutate(
    dia = factor(dia, levels = c("MCI", "AD"))
  ) |>
  ggplot(aes(
    x = dep_med,
    y = nmut,
    color = dia
  )) +
  geom_point() +
  geom_point(aes(color = dia), show.legend = FALSE) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linetype = 21) +
  geom_smooth(aes(color = dia), method = "glm", se = FALSE) +
  annotate(
    geom = "text",
    x = 500,
    y = 30,
    size = 6,
    label = latex2exp::TeX(glue::glue("$\\rho$={round(ct$estimate, 2)}, $P$={round(ct$p.value,3)}")),
    fontface = "bold"
  ) +
  ggsci::scale_color_jama(
    name = "Disease type"
  ) +
  theme_bw() +
  theme(
    # panel.grid = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 16, face = "bold", colour = "black"),
    legend.position = "bottom",
    plot.title = element_text(
      hjust = 0.5,
      size = 16,
      face = "bold",
      color = "black"
    )
  ) +
  labs(
    x = "Median depth",
    y = "# of variants",
    title = "Olah et al, Nat Commun, 2020"
  ) ->
p_nmut_median_depth
p_nmut_median_depth

ggsave(
  filename = "All-factor-correlations-linear-depth-nvariant.pdf",
  plo = p_nmut_median_depth,
  device = "pdf",
  width = 9,
  height = 6,
  path = "/home/liuc9/github/scMOCHA/03-ADKP/output"
)


# Age mutation ------------------------------------------------------------

cor.test(
  formula = ~ nmut + Age,
  data = metadata_anno_depth_dep
) ->
cta

cor.test(
  formula = ~ nmut + Age,
  data = metadata_anno_depth_dep |>
    dplyr::filter(dia == "MCI")
) ->
cta_mci

cor.test(
  formula = ~ nmut + Age,
  data = metadata_anno_depth_dep |>
    dplyr::filter(dia == "AD")
) ->
cta_ad

yhight <- 32
xwidth <- 78

metadata_anno_depth_dep |>
  # dplyr::filter(nmut >10) |>
  dplyr::mutate(
    label = glue::glue(
      "N variants = {nmut}\n Median depth = {dep_med}\n Gender = {Sex}"
    )
  ) |>
  dplyr::mutate(
    dia = factor(dia, levels = c("MCI", "AD"))
  ) |>
  ggplot(aes(
    x = Age,
    y = nmut
  )) +
  geom_point(aes(color = dia), show.legend = FALSE) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linetype = 21) +
  geom_smooth(aes(color = dia), method = "glm", se = FALSE) +
  ggrepel::geom_text_repel(
    aes(label = label),
    # box.padding = 0.5,
    max.overlaps = 10,
    # max.overlaps = Inf
    size = 3,
    min.segment.length = 0,
    seed = 42,
    box.padding = 0.5
  ) +
  ggsci::scale_color_jama(
    name = "Disease type"
  ) +
  annotate(
    geom = "segment",
    x = xwidth,
    y = yhight,
    xend = xwidth + 1,
    yend = yhight,
    linetype = 21,
    colour = "black",
    linewidth = 1
  ) +
  annotate(
    geom = "text",
    x = xwidth + 3,
    y = yhight,
    size = 5,
    label = latex2exp::TeX(glue::glue("$\\rho$={round(cta$estimate, 2)}, $P$={round(cta$p.value,3)}")),
    fontface = "bold",
  ) +
  annotate(
    geom = "segment",
    x = xwidth,
    y = yhight - 2,
    xend = xwidth + 1,
    yend = yhight - 2,
    linetype = 1,
    colour = ggsci::pal_jama()(2)[1],
    linewidth = 1
  ) +
  annotate(
    geom = "text",
    x = xwidth + 3,
    y = yhight - 2,
    size = 5,
    label = latex2exp::TeX(glue::glue("$\\rho$={round(cta_mci$estimate, 2)}, $P$={round(cta_mci$p.value,3)}")),
    fontface = "bold",
    color = ggsci::pal_jama()(2)[1],
  ) +
  annotate(
    geom = "segment",
    x = xwidth,
    y = yhight - 4,
    xend = xwidth + 1,
    yend = yhight - 4,
    linetype = 1,
    colour = ggsci::pal_jama()(2)[2],
    linewidth = 1
  ) +
  annotate(
    geom = "text",
    x = xwidth + 3,
    y = yhight - 4,
    size = 5,
    label = latex2exp::TeX(glue::glue("$\\rho$={round(cta_ad$estimate, 2)}, $P$={round(cta_ad$p.value,3)}")),
    fontface = "bold",
    color = ggsci::pal_jama()(2)[2]
  ) +
  theme_bw() +
  theme(
    # panel.grid = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 16, face = "bold", colour = "black"),
    legend.position = "top"
  ) +
  labs(
    x = "Age",
    y = "# of variants"
  ) ->
p_linear_1
p_linear_1

ggsave(
  filename = "All-factor-correlations-linear-age-nvariant.pdf",
  plo = p_linear_1,
  device = "pdf",
  width = 9,
  height = 6,
  path = "/home/liuc9/github/scMOCHA/03-ADKP/output"
)



cor.test(
  formula = ~ nmut + Age,
  data = metadata_anno_depth_dep |>
    dplyr::filter(dep_med > 1000)
) ->
cta

cor.test(
  formula = ~ nmut + Age,
  data = metadata_anno_depth_dep |>
    dplyr::filter(dep_med > 1000) |>
    dplyr::filter(dia == "MCI")
) ->
cta_mci

cor.test(
  formula = ~ nmut + Age,
  data = metadata_anno_depth_dep |>
    dplyr::filter(dep_med > 1000) |>
    dplyr::filter(dia == "AD")
) ->
cta_ad

yhight <- 32
xwidth <- 78


metadata_anno_depth_dep |>
  dplyr::filter(dep_med > 1000) |>
  dplyr::mutate(
    label = glue::glue(
      "N variants = {nmut}\n Median depth = {dep_med}\n Gender = {Sex}"
    )
  ) |>
  dplyr::mutate(
    dia = factor(dia, levels = c("MCI", "AD"))
  ) |>
  ggplot(aes(
    x = Age,
    y = nmut
  )) +
  geom_point(aes(color = dia), show.legend = FALSE) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linetype = 21) +
  geom_smooth(aes(color = dia), method = "glm", se = FALSE) +
  ggrepel::geom_text_repel(
    aes(label = label),
    # box.padding = 0.5,
    max.overlaps = 10,
    # max.overlaps = Inf
    size = 3,
    min.segment.length = 0,
    seed = 42,
    box.padding = 0.5
  ) +
  ggsci::scale_color_jama(
    name = "Disease type"
  ) +
  annotate(
    geom = "segment",
    x = xwidth,
    y = yhight,
    xend = xwidth + 1,
    yend = yhight,
    linetype = 21,
    colour = "black",
    linewidth = 1
  ) +
  annotate(
    geom = "text",
    x = xwidth + 3,
    y = yhight,
    size = 5,
    label = latex2exp::TeX(glue::glue("$\\rho$={round(cta$estimate, 2)}, $P$={round(cta$p.value,3)}")),
    fontface = "bold"
  ) +
  annotate(
    geom = "segment",
    x = xwidth,
    y = yhight - 2,
    xend = xwidth + 1,
    yend = yhight - 2,
    linetype = 1,
    colour = ggsci::pal_jama()(2)[1],
    linewidth = 1
  ) +
  annotate(
    geom = "text",
    x = xwidth + 3,
    y = yhight - 2,
    size = 5,
    label = latex2exp::TeX(glue::glue("$\\rho$={round(cta_mci$estimate, 2)}, $P$={round(cta_mci$p.value,3)}")),
    fontface = "bold",
    color = ggsci::pal_jama()(2)[1]
  ) +
  annotate(
    geom = "segment",
    x = xwidth,
    y = yhight - 4,
    xend = xwidth + 1,
    yend = yhight - 4,
    linetype = 1,
    colour = ggsci::pal_jama()(2)[2],
    linewidth = 1
  ) +
  annotate(
    geom = "text",
    x = xwidth + 3,
    y = yhight - 4,
    size = 5,
    label = latex2exp::TeX(glue::glue("$\\rho$={round(cta_ad$estimate, 2)}, $P$={round(cta_ad$p.value,3)}")),
    fontface = "bold",
    color = ggsci::pal_jama()(2)[2]
  ) +
  theme_bw() +
  theme(
    # panel.grid = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 16, face = "bold", colour = "black"),
    legend.position = "top"
  ) +
  labs(
    x = "Age",
    y = "# of variants"
  ) ->
p_linear_2
p_linear_2

ggsave(
  filename = "All-factor-correlations-linear-age-nvariant2.pdf",
  plo = p_linear_2,
  device = "pdf",
  width = 9,
  height = 6,
  path = "/home/liuc9/github/scMOCHA/03-ADKP/output"
)



t.test(nmut ~ dia, data = metadata_anno_depth_dep) |> report::report()
t.test(nmut ~ Sex, data = metadata_anno_depth_dep) |> report::report()

# footer ------------------------------------------------------------------
# future::plan(future::sequential)

# save image --------------------------------------------------------------
# save.image(file = "/scr1/users/liuc9/mitochondrial/realdata/03-ADKP/output/rda/07-integrated-analysis.rda")

# load(file = "/scr1/users/liuc9/mitochondrial/realdata/03-ADKP/output/rda/07-integrated-analysis.rda")
