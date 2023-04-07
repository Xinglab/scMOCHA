# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Mar 30 18:05:33 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(Seurat)
library(Azimuth)


# args --------------------------------------------------------------------

args <- commandArgs(TRUE)

h5file <- args[1]
refname <- args[2]
celllevel <- args[3]
# h5file <- "/scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/flu2/Flu2/outs/filtered_feature_bc_matrix.h5"
# refname <- "pbmcref"
# celllevel <- "celltype.l1"


# src ---------------------------------------------------------------------

pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv")


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------

fn_load_sc_10x <- function(.x, .project = "singlecell") {
  # .x is the feature_bc_matrix
  # .x <- feature_path
  # .x <- h5file
  # .x
  
  .counts <- tryCatch(
    expr = {
      Seurat::Read10X_h5(filename = .x)
    },
    error = function(e){
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
  
  apply(
    .sc@assays$RNA@counts,
    2,
    function(x) (100 * max(x)) / sum(x)
  ) ->
    .sc$Percent.Largest.Gene
  
  
  .sc@meta.data %>%
    dplyr::arrange(percent.mt) %>%
    ggplot(aes(nCount_RNA, nFeature_RNA, color = percent.mt)) +
    geom_point() +
    scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
    ggtitle("Mito of plotting QC metrics") +
    geom_hline(yintercept = 500, color = "red") +
    geom_hline(yintercept = 6000, color = "red") +
    theme_bw() ->
    .metrics_mito
  
  # ggsave(
  #   plot = .metrics_mito,
  #   filename = "{.project}-metrics-mt.pdf" %>% glue::glue(),
  #   device = "pdf",
  #   width = 7,
  #   height = 5
  # )
  
  .sc@meta.data %>%
    dplyr::arrange(percent.ribo) %>%
    ggplot(aes(nCount_RNA, nFeature_RNA, color = percent.ribo)) +
    geom_point() +
    scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
    ggtitle("Ribo of plotting QC metrics") +
    geom_hline(yintercept = 500, color = "red") +
    geom_hline(yintercept = 6000, color = "red") +
    theme_bw() ->
    .metrics_ribo
  
  # ggsave(
  #   plot = .metrics_ribo,
  #   filename = "{.project}-metrics-ribo.pdf" %>% glue::glue(),
  #   device = "pdf",
  #   width = 7,
  #   height = 5
  # )
  
  .plot_2 <- (.metrics_mito + .metrics_ribo) +
    plot_annotation(
      title = glue::glue("Quality control {.project}"),
      tag_levels = "A"
    )
  
  # ggsave(
  #   plot = .plot,
  #   filename = "{.project}-metrics-mt-ribo.pdf" %>% glue::glue(),
  #   device = "pdf",
  #   width = 14,
  #   height = 5
  # )
  
  .sc@meta.data %>%
    ggplot(aes(percent.mt)) +
    geom_histogram(binwidth = 0.5, fill="red") +
    ggtitle("Distribution of Percentage Mitochondrion") +
    geom_vline(xintercept = 75) +
    theme_bw() ->
    .percent_mt
  
  .sc@meta.data %>%
    ggplot(aes(percent.ribo)) +
    geom_histogram(binwidth = 0.5, fill="green") +
    ggtitle("Distribution of Percentage Ribosome") +
    geom_vline(xintercept = 50) +
    theme_bw() ->
    .percent_ribo
  
  .sc@meta.data %>%
    ggplot(aes(Percent.Largest.Gene)) +
    geom_histogram(binwidth = 0.7, fill="blue") +
    ggtitle("Distribution of Percentage Largest Gene") +
    geom_vline(xintercept = 50) +
    theme_bw() ->
    .percent_largest_gene
  
  .plot_3 <- (.percent_mt | .percent_ribo | .percent_largest_gene) +
    plot_annotation(
      title = glue::glue("Quality control {.project}"),
      tag_levels = "A"
    )
  
  # ggsave(
  #   plot = .plot,
  #   filename = "{.project}-qc-mt-ribo-largest.pdf" %>% glue::glue(),
  #   device = "pdf",
  #   width = 15,
  #   height = 5
  # )
  
  .sc_sub <- subset(
    x = .sc,
    subset = nFeature_RNA > 500 &
      nFeature_RNA < 6000 &
      percent.mt < 75 &
      percent.ribo < 50 &
      Percent.Largest.Gene < 50
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
    `number of cells after filtering` = .n_y_cells
  )
}

fn_azimuth <- function(.sc, .ref, .celllevel) {
  # .x <- sc$sc_filter
  
  .sca <- Azimuth::RunAzimuth(
    query = .sc,
    reference = .ref
  )
  
  .celltype <-  .sca[[glue::glue("predicted.{.celllevel}")]][, 1] |> factor()

  .celltype_collapse <- gsub(
    pattern = " ",
    replacement = "_",
    x = .celltype
  ) |> factor()
  
  .sca[["celltype"]] <- .celltype
  .sca[["celltype_name"]] <- .celltype_collapse
  
  .sca
  
}

fn_plot_azimuth_umap <- function(.x) {
  
  .umap <- .x@reductions$ref.umap@cell.embeddings |> data.table::as.data.table()
  colnames(.umap) <- c("UMAP_1", "UMAP_2")
  
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
  
  .xxx_celltype |>
    dplyr::ungroup() %>% 
    dplyr::mutate(csum = rev(cumsum(rev(n)))) %>% 
    dplyr::mutate(pos = n/2 + dplyr::lead(csum, 1)) %>% 
    dplyr::mutate(pos = dplyr::if_else(is.na(pos), n/2, pos)) %>% 
    dplyr::mutate(percentage = n/sum(n)) %>% 
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
    ) ->
    .p_pie
  
  
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
      
      if(nrow(.mmm) == 1) {
        return(
          .mmd |>
            # dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) |>
            dplyr::summarise(UMAP_1  = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
        )
        
      }
      
      .fc <- .mmm$n[[1]] / .mmm$n[[2]] # 1.1
      
      if(.fc > 1.1) {
        .mmd |>
          dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) |>
          dplyr::summarise(UMAP_1  = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
      } else {
        .mmd |>
          # dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) |>
          dplyr::summarise(UMAP_1  = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
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
        override.aes = list(size=4)
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
    coord_fixed(
      ratio = 1
    ) ->
    .p
  
  
  list(
    plot_umap = .p,
    plot_celltype_pie = .p_pie,
    cell_ratio = .xxx_celltype
  )
  
}


# load data ---------------------------------------------------------------

sc <- fn_load_sc_10x(h5file)

# body --------------------------------------------------------------------


# Stat --------------------------------------------------------------------


sc$cell_stats <- fn_stat_cell(
  .x = sc$sc,
  .y = sc$sc_filter
)


# Azimuth -----------------------------------------------------------------

sc$sc_azimuth <- fn_azimuth(
  .sc = sc$sc_filter,
  .ref = refname,
  .celllevel = celllevel
)


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


sc$cellbarcode_bulk <- sc$cellbarcode_cluster |> 
  dplyr::mutate(cluster = "Bulk")


# Plot umap ---------------------------------------------------------------

sc$plot_umap <- fn_plot_azimuth_umap(.x = sc$sc_azimuth)




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

ggsave(
  filename = "plot-umap.pdf",
  plot = sc$plot_umap$plot_umap,
  device = "pdf",
  width = 9,
  height = 7
)

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
