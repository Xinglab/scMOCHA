#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Mar 30 18:05:33 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Azimuth))
library(GetoptLong)
library(logger)

# args --------------------------------------------------------------------

# args <- commandArgs(TRUE)
#
# h5file <- args[1]
# npcs <- args[2]
# reso <- args[3]
# refname <- args[4]
# celllevel <- args[5]

# h5file <- "/scr1/users/liuc9/mitochondrial/realdata/05-Liming/scmocha-mixed-cellline-high-depth/cromwell-executions/scMOCHA/023d7328-9097-4e50-8c11-19f860c5519e/call-cell_cluster_annotation/inputs/1509575042/filtered_feature_bc_matrix.h5"
# npcs <- 10
# reso <- 0.1
# refname <- "/home/liuc9/github/scMOCHA/03-ADKP/forrefs/azimuth_syn21438358"
# celllevel <- "annotation.l1"

# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.

verbose <- FALSE
npcs <- 10
reso <- 0.1
refname_celllevel <- list(
  refname = NA_character_,
  celllevel = NA_character_
)
nFeature_RNA_min <- 200
nFeature_RNA_max <- 8000
percent_mt_max <- 75
percent_ribo_max <- 50
percent_Lagest_Gene_max <- 50


spec <- "
Usage: Rscript azimuth.R [options]

Options:
<h5file=s> possorted_genome_bam.MT.depth
<npcs=f> default 10
<reso=f> default 0.1
<refname_celllevel=s%> azimuth reference name cell type level, should be used as -refname_celllevel refname=${refname} celllevel=${celllevel}
<nFeature_RNA_min=i> default 500
<nFeature_RNA_max=i> default 8000
<percent_mt_max=i> default 75
<percent_ribo_max=i> default 50
<percent_Lagest_Gene_max=i> default 50
<verbose!> Print messages
"

GetoptLong.options(help_style = "two-column")
GetoptLong(spec, template_control = list(opt_width = 50))

log_threshold(TRACE)
log_layout(layout_glue_colors)

# print(celllevel)
refname <- refname_celllevel$refname
celllevel <- refname_celllevel$celllevel

use_azimuth <- TRUE

if (!file.exists(h5file)) {
  # message(
  #   "Notice: {h5file} does not exist" |> glue::glue()
  # )
  log_error(h5file, " does not exists!")
  quit(save = "no")
} else {
  # message(
  #   "Notice: {h5file} exists" |> glue::glue()
  # )
  log_success(h5file, " exists.")
}

if (is.na(refname) || refname == "") {
  # message(
  #   "Notice: refname is not defined \n the cell cluster and annotation will not be using by Azimuth"
  # )
  log_warn("refname is not defined the cell cluster and annotation will not be used by Azimuth")
  use_azimuth <- FALSE
}

if (is.na(celllevel) || celllevel == "") {
  # message(
  #   "Notice: celllevel is not defined \n the cell cluster and annotation will not be using by Azimuth"
  # )
  log_warn("celllevel is not defined the cell cluster and annotation will not be used by Azimuth")
  use_azimuth <- FALSE
}

if (use_azimuth) {
  log_info("You will use Azimuth reference cell map by use tissue ", refname, " and the cell level is ", celllevel)
} else {
  log_warn("You will not use Azimuth reference cell map for cell annotation")
}

#
log_success("h5file ", h5file)
log_success("npcs ", npcs, " ", class(npcs))
log_success("reso ", reso, " ", class(reso))
log_success("refname ", refname, " ", class(refname))
log_success("celllevel ", celllevel, " ", class(celllevel))
log_success("nFeature_RNA_min ", nFeature_RNA_min, " ", class(nFeature_RNA_min))
log_success("nFeature_RNA_max ", nFeature_RNA_max, " ", class(nFeature_RNA_max))
log_success("percent_mt_max ", percent_mt_max, " ", class(percent_mt_max))
log_success("percent_ribo_max ", percent_ribo_max, " ", class(percent_ribo_max))
log_success("percent_Lagest_Gene_max ", percent_Lagest_Gene_max, " ", class(percent_Lagest_Gene_max))
# quit(save="no", status = 0, runLast = F)
# src ---------------------------------------------------------------------

pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv")



# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------

fn_metrics_mito <- function(.sc) {
  .sc@meta.data %>%
    dplyr::arrange(percent.mt) %>%
    ggplot(aes(nCount_RNA, nFeature_RNA, color = percent.mt)) +
    geom_point() +
    scale_color_gradientn(colors = c("black", "blue", "green2", "red", "yellow")) +
    ggtitle("Mito of plotting QC metrics") +
    geom_hline(yintercept = 500, color = "red") +
    geom_hline(yintercept = 6000, color = "red") +
    theme_bw() ->
  .metrics_mito
}

fn_metrics_ribo <- function(.sc) {
  .sc@meta.data %>%
    dplyr::arrange(percent.ribo) %>%
    ggplot(aes(nCount_RNA, nFeature_RNA, color = percent.ribo)) +
    geom_point() +
    scale_color_gradientn(colors = c("black", "blue", "green2", "red", "yellow")) +
    ggtitle("Ribo of plotting QC metrics") +
    geom_hline(yintercept = 500, color = "red") +
    geom_hline(yintercept = 6000, color = "red") +
    theme_bw() ->
  .metrics_ribo
}

fn_percent_mt_ribo_lg <- function(.sc) {
  .sc@meta.data %>%
    ggplot(aes(percent.mt)) +
    geom_histogram(binwidth = 0.5, fill = "red") +
    ggtitle("Distribution of Percentage Mitochondrion") +
    geom_vline(xintercept = 75) +
    theme_bw() ->
  .percent_mt

  .sc@meta.data %>%
    ggplot(aes(percent.ribo)) +
    geom_histogram(binwidth = 0.5, fill = "green") +
    ggtitle("Distribution of Percentage Ribosome") +
    geom_vline(xintercept = 50) +
    theme_bw() ->
  .percent_ribo

  .sc@meta.data %>%
    ggplot(aes(Percent.Largest.Gene)) +
    geom_histogram(binwidth = 0.7, fill = "blue") +
    ggtitle("Distribution of Percentage Largest Gene") +
    geom_vline(xintercept = 50) +
    theme_bw() ->
  .percent_largest_gene

  (.percent_mt | .percent_ribo | .percent_largest_gene)
}

fn_create_sc <- function(.x, .project = "singlecell") {
  .counts <- tryCatch(
    expr = {
      Seurat::Read10X_h5(filename = .x)
    },
    error = function(e) {
      .xx <- gsub(pattern = ".h5", replacement = "", x = .x)
      Seurat::Read10X(data.dir = .xx)
    }
  )
  # .counts <- Seurat::Read10X_h5(filename = .x)

  .sc <- Seurat::CreateSeuratObject(
    counts = .counts,
    project = .project,
    # min.cells = 3,
    # min.features = 200
  )

  .sc <- Seurat::PercentageFeatureSet(
    object = .sc,
    pattern = "^MT-",
    col.name = "percent.mt"
  )

  .sc <- Seurat::PercentageFeatureSet(
    object = .sc,
    pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",
    # pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa",
    col.name = "percent.ribo"
  )

  if (packageVersion("Seurat") > "5") {
    apply(
      # .sc@assays$RNA@counts, # Seurat version 4
      .sc@assays$RNA@layers$counts, # Seurat version 5
      2,
      function(x) (100 * max(x)) / sum(x)
    ) ->
    .sc$Percent.Largest.Gene
  } else {
    apply(
      .sc@assays$RNA@counts, # Seurat version 4
      # .sc@assays$RNA@layers$counts, # Seurat version 5
      2,
      function(x) (100 * max(x)) / sum(x)
    ) ->
    .sc$Percent.Largest.Gene
  }

  .sc
}

fn_load_sc_10x <- function(.x, .project = "singlecell") {
  .sc <- fn_create_sc(.x, .project)

  .plot_2 <- (fn_metrics_mito(.sc) + fn_metrics_ribo(.sc)) +
    plot_annotation(
      title = glue::glue("Quality control {.project}"),
      tag_levels = "A"
    )

  .plot_3 <- fn_percent_mt_ribo_lg(.sc) +
    plot_annotation(
      title = glue::glue("Quality control {.project}"),
      tag_levels = "A"
    )

  .sc_sub <- subset(
    x = .sc,
    subset = nFeature_RNA > nFeature_RNA_min &
      nFeature_RNA < nFeature_RNA_max &
      percent.mt < percent_mt_max &
      percent.ribo < percent_ribo_max &
      Percent.Largest.Gene < percent_Lagest_Gene_max
  )

  list(
    sc = .sc,
    sc_filter = .sc_sub,
    plot_metrics = .plot_2,
    plot_qc = .plot_3
  )
}

fn_stat_cell <- function(.x, .y) {
  # .x <- project_sct$sc[[1]]
  # .y <- project_sct$sct[[1]]

  .xd <- .x@meta.data
  .yd <- .y@meta.data

  .n_x_cells <- nrow(.xd)
  .median_umicount_per_cell <- median(.xd$nCount_RNA)
  .median_gene_per_cell <- median(.xd$nFeature_RNA)
  .n_y_cells <- nrow(.yd)

  tibble::tibble(
    `estimated number of cells` = .n_x_cells,
    # `mean reads per cell` = .mean_reads_per_cell,
    `median UMI counts per cell` = .median_umicount_per_cell,
    `median genes per cell` = .median_gene_per_cell,
    `number of cells after filtering` = .n_y_cells,
    `Preserve ratio` = .n_y_cells / .n_x_cells
  )
}

fn_azimuth <- function(.sc, .ref, .celllevel) {
  # .sc <- sc$sc_filter
  # .sc <- sc$sc

  .sca <- Azimuth::RunAzimuth(
    query = .sc,
    reference = .ref
  )

  .celltype <- .sca[[glue::glue("predicted.{.celllevel}")]][, 1] |> factor()

  .celltype_collapse <- gsub(
    pattern = "[[:punct:]]| ",
    replacement = "_",
    x = .celltype
  ) |> factor()

  .sca[["celltype"]] <- .celltype
  .sca[["celltype_name"]] <- .celltype_collapse

  .sca
}

fn_scnorm <- function(.sc) {
  .sc |> Seurat::NormalizeData() -> .scn

  .scn <- Seurat::FindVariableFeatures(
    .scn,
    selection.method = "vst",
    nfeatures = 2000
  )

  .allgenes <- rownames(.scn)
  .scn <- Seurat::ScaleData(
    .scn,
    features = .allgenes
  )

  .npcs <- as.numeric(npcs)
  .reso <- as.numeric(reso)

  .scn |>
    Seurat::RunPCA(features = VariableFeatures(.scn)) |>
    Seurat::FindNeighbors(reduction = "pca", dims = .npcs) |>
    Seurat::FindClusters(resolution = .reso) |>
    Seurat::RunUMAP(reduction = "pca", dims = 1:.npcs) |>
    Seurat::RunTSNE(reduction = "pca", dims = 1:.npcs) ->
  .scna

  .celltype <- glue::glue("cluster_{.scna[['seurat_clusters']][, 1]}") |> factor()
  .celltype_collapse <- .celltype

  .scna[["celltype"]] <- .celltype
  .scna[["celltype_name"]] <- .celltype_collapse

  .scna
}

fn_sctransform <- function(.sc) {
  .sct <- Seurat::SCTransform(
    object = .sc,
    vars.to.regress = c("percent.mt", "percent.ribo")
  )

  .npcs <- as.numeric(npcs)
  .reso <- as.numeric(reso)

  .sct |>
    Seurat::RunPCA(dim = .npcs) |>
    Seurat::FindNeighbors(reduction = "pca", dims = 1:.npcs) |>
    Seurat::FindClusters(resolution = 0.1) |>
    Seurat::RunUMAP(reduction = "pca", dims = 1:(.npcs / 2)) |>
    Seurat::RunTSNE(reduction = "pca", dims = 1:.npcs) ->
  .scta

  .celltype <- glue::glue("cluster_{.scta[['seurat_clusters']][, 1]}") |> factor()
  .celltype_collapse <- .celltype

  .scta[["celltype"]] <- .celltype
  .scta[["celltype_name"]] <- .celltype_collapse

  # fn_plot_azimuth_umap(.scta, .plottype = "tsne")$plot_umap

  .scta
}

fn_cluster_anno <- function(.sc, .use_azimuth, .ref, .celllevel) {
  # .sc <- sc$sc_filter
  # .sc <- sc$sc

  if (.use_azimuth) {
    .sca <-
      tryCatch(
        expr = {
          fn_azimuth(.sc, .ref, .celllevel)
        },
        error = \(e) {
          fn_sctransform(.sc)
        }
      )
  } else {
    .sca <- fn_sctransform(.sc)
  }

  .sca
}

fn_plot_cluster <- function(.xxx, .xy_labels = c("UMAP_1", "UMAP_2")) {
  .xxx |>
    dplyr::group_by(celltype) |>
    tidyr::nest() |>
    dplyr::ungroup() |>
    dplyr::mutate(u = purrr::map(.x = data, .f = function(.m) {
      # d |>
      #   dplyr::filter(cluster == 14) |>
      #   dplyr::pull(data) |>
      #   .[[1]] ->
      #   .m

      .m |>
        dplyr::summarise(u1 = mean(UMAP_1), u2 = mean(UMAP_2)) ->
      .mm

      .m |>
        dplyr::mutate(u1 = UMAP_1 > .mm$u1, u2 = UMAP_2 > .mm$u2) ->
      .mmd

      .mmd |>
        dplyr::group_by(u1, u2) |>
        dplyr::count() |>
        dplyr::ungroup() |>
        dplyr::arrange(-n) ->
      .mmm

      if (nrow(.mmm) == 1) {
        return(
          .mmd |>
            # dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) |>
            dplyr::summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
        )
      }

      .fc <- .mmm$n[[1]] / .mmm$n[[2]] # 1.1

      if (.fc > 1.1) {
        .mmd |>
          dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) |>
          dplyr::summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
      } else {
        .mmd |>
          # dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) |>
          dplyr::summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
      }
    })) |>
    dplyr::select(-data) |>
    tidyr::unnest(cols = u) |>
    dplyr::arrange(celltype) ->
  .xxx_label

  ggplot() +
    geom_point(
      data = .xxx,
      aes(
        x = UMAP_1,
        y = UMAP_2,
        colour = celltype,
        shape = NULL,
        alpha = NULL
      ),
      size = 0.7
    ) +
    geom_text(
      data = .xxx_label,
      aes(
        label = celltype,
        x = UMAP_1,
        y = UMAP_2,
      ),
      size = 6
    ) +
    scale_colour_manual(
      name = "Cell type",
      values = pcc$color,
      labels = .xxx_label$celltype,
      guide = guide_legend(
        ncol = 1,
        override.aes = list(size = 4)
      )
    ) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(
        colour = "black",
        linewidth = 0.5,
        arrow = grid::arrow(
          angle = 5,
          length = unit(5, "npc"),
          type = "closed"
        )
      ),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.05
      ),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(
        face = "bold",
        color = "black",
        size = 14
      ),
      legend.text = element_text(
        face = "bold",
        color = "black",
        size = 12
      )
    ) +
    labs(
      x = .xy_labels[[1]],
      y = .xy_labels[[2]]
    ) +
    coord_fixed(
      ratio = 1
    )
}

fn_celltype_pie_plot <- function(.xxx_celltype) {
  .xxx_celltype |>
    dplyr::ungroup() %>%
    dplyr::mutate(csum = rev(cumsum(rev(n)))) %>%
    dplyr::mutate(pos = n / 2 + dplyr::lead(csum, 1)) %>%
    dplyr::mutate(pos = dplyr::if_else(is.na(pos), n / 2, pos)) %>%
    dplyr::mutate(percentage = n / sum(n)) %>%
    ggplot(aes(x = "", y = n, fill = celltype)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    # scale_fill_brewer(palette = "Dark2", name = NULL) +
    scale_fill_manual(
      name = NULL,
      values = pcc$color,
    ) +
    # scale_color_manual(
    #   name = NULL,
    #   values = pcc$color
    # ) +
    ggrepel::geom_label_repel(
      aes(
        y = pos,
        label = glue::glue("{celltype}\n{n} ({scales::percent(percentage)})"),
        fill = celltype,
        # color = celltype
      ),
      size = 6,
      # fill = "white",
      nudge_x = 1,
      show.legend = FALSE,
    ) +
    coord_polar(theta = "y", start = 0) +
    theme_void() +
    theme(
      plot.title = element_text(
        # vjust = -2,
        hjust = 0.5,
        size = 22,
      ),
      legend.position = "none"
    )
}

fn_plot_azimuth_umap <- function(.x) {
  .col_names <- c("UMAP_1", "UMAP_2")

  if ("ref.umap" %in% names(.x@reductions)) {
    .umap <- .x@reductions$ref.umap@cell.embeddings |> data.table::as.data.table()
    colnames(.umap) <- .col_names
    .tsne <- NULL
  } else {
    .umap <- .x@reductions$umap@cell.embeddings |> data.table::as.data.table()
    colnames(.umap) <- .col_names
    .tsne <- .x@reductions$tsne@cell.embeddings |> data.table::as.data.table()
    colnames(.tsne) <- .col_names
  }

  # .umap
  .x@meta.data |>
    dplyr::select(
      celltype
    ) |>
    data.table::as.data.table() ->
  .xx

  .xxx <- dplyr::bind_cols(.umap, .xx)

  .xxx |>
    dplyr::group_by(celltype) |>
    dplyr::count() |>
    dplyr::ungroup() |>
    dplyr::mutate(ratio = n / sum(n)) ->
  .xxx_celltype

  .plot_pie <- fn_celltype_pie_plot(.xxx_celltype)

  .plot_umap <- fn_plot_cluster(.xxx)
  .plot_tsne <- if (is.null(.tsne)) {
    NULL
  } else {
    .xxx <- dplyr::bind_cols(.tsne, .xx)
    fn_plot_cluster(.xxx, .xy_labels = c("tSNE_1", "tSNE_2"))
  }

  list(
    plot_umap = .plot_umap,
    plot_tsne = .plot_tsne,
    plot_celltype_pie = .plot_pie,
    cell_ratio = .xxx_celltype
  )
}

fn_check_cellref <- function(.refname) {
  # SeuratData::InstalledData() |> dplyr::glimpse()
  if (dir.exists(.refname)) {
    message(glue::glue("Azimuth reference {.refname} installed"))
    use_azimuth <<- TRUE
    return(1)
  }

  .sd <- SeuratData::AvailableData() |>
    dplyr::filter(
      grepl("Azimuth Reference", x = Summary)
    )

  .ref <- .sd |>
    dplyr::filter(Dataset == .refname)

  if (length(.ref$Installed) != 0 && .ref$Installed) {
    message(glue::glue("Azimuth reference {.refname} installed"))
  } else {
    tryCatch(
      expr = {
        SeuratData::InstallData(
          ds = .refname
        )
      },
      warning = \(w) {
        use_azimuth <<- FALSE
      },
      error = \(e) {
        use_azimuth <<- FALSE
      }
    )
  }
}


# load data ---------------------------------------------------------------
log_warn("Check refname exists")
fn_check_cellref(refname)


# Load 10x ----------------------------------------------------------------

sc <- fn_load_sc_10x(h5file)
log_success(h5file, " loaded!!!")

# body --------------------------------------------------------------------


# Stat --------------------------------------------------------------------


sc$cell_stats <- fn_stat_cell(
  .x = sc$sc,
  .y = sc$sc_filter
)
log_success("Stats done!!!")


# Azimuth -----------------------------------------------------------------

# message("Notice: fn_cluster_anno is running")
log_info("Notice: fn_cluster_anno is running")

sc$sc_azimuth <- fn_cluster_anno(
  .sc = sc$sc_filter,
  .use_azimuth = use_azimuth,
  .ref = refname,
  .celllevel = celllevel
)
log_success("Notice: fn_cluster_anno is done")



# Cell barcode ------------------------------------------------------------

sc$sc_azimuth@meta.data |>
  tibble::rownames_to_column(
    var = "cellbarcode"
  ) |>
  data.table::as.data.table() |>
  dplyr::select(
    cellbarcode,
    celltype_name
  ) |>
  dplyr::mutate(
    tag = "CJ",
    cluster = celltype_name
  ) |>
  dplyr::select(1, 3, 4) ->
sc$cellbarcode_cluster
log_success("barcode")
# quit(save = "no")

sc$cellbarcode_bulk <- sc$cellbarcode_cluster |>
  dplyr::mutate(cluster = "Bulk")

log_success("barcode bulk")
# Plot umap ---------------------------------------------------------------

sc$plot_umap <- fn_plot_azimuth_umap(.x = sc$sc_azimuth)

log_success("Plot umap!")


# readr -------------------------------------------------------------------

# names(sc)

readr::write_tsv(
  x = sc$cellbarcode_cluster,
  file = "barcode_cluster.tsv",
  col_names = F
)
readr::write_tsv(
  x = sc$cellbarcode_bulk,
  file = "barcode_bulk.tsv",
  col_names = F
)
readr::write_tsv(
  x = sc$plot_umap$cell_ratio,
  file = "celltype_ratio.tsv"
)
writexl::write_xlsx(
  x = sc$cell_stats,
  path = "qc_cell_stats.xlsx"
)



# save plot ---------------------------------------------------------------


names(sc)
ggsave(
  filename = "plot-metrics.pdf",
  plot = sc$plot_metrics,
  device = "pdf",
  width = 12,
  height = 5
)

ggsave(
  filename = "plot-qc.pdf",
  plot = sc$plot_qc,
  device = "pdf",
  width = 12,
  height = 5
)

if (is.null(sc$plot_umap$plot_tsne)) {
  ggsave(
    filename = "plot-umap.pdf",
    plot = sc$plot_umap$plot_umap,
    device = "pdf",
    width = 9,
    height = 7
  )
} else {
  ggsave(
    filename = "plot-umap.pdf",
    plot = sc$plot_umap$plot_umap / sc$plot_umap$plot_tsne,
    device = "pdf",
    width = 9,
    height = 12
  )
}


ggsave(
  filename = "plot-pie-celltype.pdf",
  plot = sc$plot_umap$plot_celltype_pie,
  device = "pdf",
  width = 8,
  height = 5
)


# save --------------------------------------------------------------------

readr::write_rds(
  x = sc,
  file = "sc_azimuth.rds.gz",
  compress = "none"
)

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------

save.image(file = "azimuth.rda")
