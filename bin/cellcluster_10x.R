# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Nov 29 02:29:13 2022
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(Seurat)

# args --------------------------------------------------------------------

args <- commandArgs(TRUE)

h5file <- args[1]
# h5file <- "/scr1/users/liuc9/tmp/singlecell/pbmc_10k_v3/outs/filtered_feature_bc_matrix.h5"
h5file <- "/scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/flu2/Flu2/outs/raw_feature_bc_matrix.h5"


# src ---------------------------------------------------------------------
pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv")

# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

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
    min.cells = 3,
    min.features = 200
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

  ggsave(
    plot = .metrics_mito,
    filename = "{.project}-metrics-mt.pdf" %>% glue::glue(),
    device = "pdf",
    width = 7,
    height = 5
  )

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

  ggsave(
    plot = .metrics_ribo,
    filename = "{.project}-metrics-ribo.pdf" %>% glue::glue(),
    device = "pdf",
    width = 7,
    height = 5
  )

  .plot <- (.metrics_mito + .metrics_ribo) +
    plot_annotation(
      title = glue::glue("Quality control {.project}"),
      tag_levels = "A"
    )
  ggsave(
    plot = .plot,
    filename = "{.project}-metrics-mt-ribo.pdf" %>% glue::glue(),
    device = "pdf",
    width = 14,
    height = 5
  )

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

  .plot <- (.percent_mt | .percent_ribo | .percent_largest_gene) +
    plot_annotation(
      title = glue::glue("Quality control {.project}"),
      tag_levels = "A"
    )

  ggsave(
    plot = .plot,
    filename = "{.project}-qc-mt-ribo-largest.pdf" %>% glue::glue(),
    device = "pdf",
    width = 15,
    height = 5
  )

  .sc
}

fn_filter_sct <- function(.sc) {

  .sc_sub <- subset(
    x = .sc,
    subset = nFeature_RNA > 500 &
      nFeature_RNA < 6000 &
      percent.mt < 75 &
      percent.ribo < 50 &
      Percent.Largest.Gene < 50
  )

  .sc_sub_sct <- Seurat::SCTransform(
    object = .sc_sub,
    do.scale = FALSE,
    do.center = FALSE
  )
  .sc_sub_sct <- Seurat::CellCycleScoring(
    object = .sc_sub_sct,
    g2m.features = Seurat::cc.genes$g2m.genes,
    s.features = Seurat::cc.genes$s.genes
  )
  .sc_sub_sct$CC.Difference <- .sc_sub_sct$S.Score - .sc_sub_sct$G2M.Score

  Seurat::DefaultAssay(.sc_sub_sct) <- "RNA"

  .sc_sub_sct_sct <- Seurat::SCTransform(
    object = .sc_sub_sct,
    method = "glmGamPoi",
    vars.to.regress = c("percent.mt", "percent.ribo", "CC.Difference"),
    do.scale = TRUE,
    do.center = TRUE
  )

  .sc_sub_sct_sct

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

fn_plot_umap <- function(.x, .celltype="sctype", .reduction="umap") {
  .x <- sct_cluster
  # .celltype="sctype"
  # .reduction="umap"
  #
  .umap <- as.data.frame(.x@reductions[[.reduction]]@cell.embeddings)
  colnames(.umap) <- c("UMAP_1", "UMAP_2")
  .xx <- .x@meta.data[, c("seurat_clusters", .celltype)] %>%
    dplyr::rename(cluster = seurat_clusters, celltype =  .celltype)
  .xxx <- dplyr::bind_cols(.umap, .xx)
  .xxx %>%
    dplyr::select(cluster, celltype) %>%
    dplyr::group_by(cluster, celltype) %>%
    dplyr::count() %>%
    dplyr::arrange(-n) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-n) %>%
    dplyr::mutate(celltype = glue::glue("{cluster} {celltype}")) ->
    .xxx_celltype

  .xxx %>%
    dplyr::group_by(cluster) %>%
    tidyr::nest() %>%
    dplyr::mutate(u = purrr::map(.x = data, .f = function(.m) {
      # d %>%
      #   dplyr::filter(cluster == 14) %>%
      #   dplyr::pull(data) %>%
      #   .[[1]] ->
      #   .m

      .m %>%
        dplyr::summarise(u1 = mean(UMAP_1), u2 = mean(UMAP_2)) ->
        .mm

      .m %>%
        dplyr::mutate(u1 = UMAP_1 > .mm$u1, u2 = UMAP_2 > .mm$u2) ->
        .mmd

      .mmd %>%
        dplyr::group_by(u1, u2) %>%
        dplyr::count() %>%
        dplyr::ungroup() %>%
        dplyr::arrange(-n) ->
        .mmm

      .fc <- .mmm$n[[1]] / .mmm$n[[2]] # 1.1
      .mmm
      .fc

      if(.fc > 1.1) {
        .mmd %>%
          dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) %>%
          dplyr::summarise(UMAP_1  = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
      } else {
        .mmd %>%
          # dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) %>%
          dplyr::summarise(UMAP_1  = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
      }

    })) %>%
    dplyr::ungroup() %>%
    tidyr::unnest(cols = u) %>%
    dplyr::select(-data) %>%
    dplyr::left_join(.xxx_celltype, by = "cluster") %>%
    dplyr::arrange(cluster) ->
    .xxx_label

  # .xxx %>%
  #   dplyr::group_by(cluster) %>%
  #   dplyr::summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2)) %>%
  #   dplyr::left_join(.xxx_celltype, by = "cluster")
  #


  .labs <- if (.reduction == "umap") {
    labs(
      x = "UMAP1",
      y = "UMAP2"
    )
  } else {
    labs(
      x = "tSNE1",
      y = "tSNE2"
    )
  }

  ggplot() +
    geom_point(
      data = .xxx,
      aes(
        x = UMAP_1,
        y = UMAP_2,
        colour = cluster,
        shape = NULL,
        alpha = NULL
      ),
      size = 0.7
    ) +
    geom_text(
      data = .xxx_label,
      aes(
        label = cluster,
        x = UMAP_1,
        y = UMAP_2,
      ),
      size = 6
    ) +
    scale_colour_manual(
      name = NULL,
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
        size = 0.5,
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
      legend.text = element_text(
        face = "bold",
        color = "black",
        size = 10
      )
    ) +
    coord_fixed(
      ratio = 1,
    ) +
    .labs
}

fn_gene_dotplot <- function(.sct_cluster, .marker, .n = 3) {

  .marker %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(n = .n, order_by = avg_log2FC) %>%
    print(n = Inf) ->
    .marker_head

  DefaultAssay(.sct_cluster) <- "SCT"

  DotPlot(
    .sct_cluster,
    features = unique(.marker_head$gene),
    cols = c("blue", "red"),
    dot.scale = 8
  ) +
    RotatedAxis()
}

# load data ---------------------------------------------------------------

# body --------------------------------------------------------------------


# sct ---------------------------------------------------------------------



sc <- fn_load_sc_10x(h5file)
readr::write_rds(
  x = sc,
  file = "sc_cluster.rds.gz"
)

sct <- fn_filter_sct(.sc = sc)


cell_stats <- fn_stat_cell(sc, sct)
writexl::write_xlsx(
  x = cell_stats,
  path = "reads-stat.xlsx"
)


# pca ---------------------------------------------------------------------


sct %>%
  Seurat::RunPCA(npcs = 30) %>%
  Seurat::RunUMAP(reduction = "pca", dims = 1:10) %>%
  Seurat::RunTSNE(reduction = "pca", dims = 1:10) %>%
  Seurat::FindNeighbors(reduction = "pca", dims = 1:10) %>%
  Seurat::FindClusters(resolution = 0.2) ->
  sct_cluster

readr::write_rds(
  x = sct_cluster,
  file = "sct_cluster.rds.gz"
)


# cell cluster ------------------------------------------------------------

sct_cluster@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column(
    var = "cellbarcode"
  ) %>%
  dplyr::select(
    cellbarcode,
    seurat_clusters
  ) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(tag = "CJ") %>%
  dplyr::mutate(seurat_clusters = as.numeric(seurat_clusters)) %>%
  dplyr::mutate(
    cluster = purrr::map_chr(
      .x = seurat_clusters,
      .f = function(.x) {
        glue::glue("cluster{.x}-1")
      }
    )
  ) %>%
  dplyr::select(1, 3, 4) ->
  cellbarcode


cellbarcode %>%
  readr::write_tsv(
    file = "barcode_cluster.tsv",
    col_names = F
  )

cellbarcode %>%
  dplyr::mutate(cluster = "Bulk") %>%
  readr::write_tsv(
    file = "barcode_bulk.tsv",
    col_names = F
  )

# annotation --------------------------------------------------------------

# sc-Type -----------------------------------------------------------------
library(HGNChelper)
# load gene set preparation function
source("https://raw.githubusercontent.com/chunjie-sam-liu/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/chunjie-sam-liu/sc-type/master/R/sctype_score_.R")
# DB file
db_full = "https://raw.githubusercontent.com/chunjie-sam-liu/sc-type/master/ScTypeDB_full.xlsx";
db_short <- "https://raw.githubusercontent.com/chunjie-sam-liu/sc-type/master/ScTypeDB_short.xlsx"
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
gs_list_immune_system = gene_sets_prepare(db_full, "Immune system")

es.max <- sctype_score(
  scRNAseqData = sct_cluster[["SCT"]]@scale.data,
  scaled = TRUE,
  gs = gs_list_immune_system$gs_positive,
  gs2 = gs_list_immune_system$gs_negative
)

cL_results <- do.call(
  "rbind",
  lapply(
    unique(sct_cluster@meta.data$seurat_clusters),
    function(cl) {
      es.max.cl <- sort(
        rowSums(es.max[, rownames(sct_cluster@meta.data[sct_cluster@meta.data$seurat_clusters == cl, ])]),
        decreasing = TRUE
      )
      head(
        data.frame(
          cluster = cl,
          type = names(es.max.cl),
          scores = es.max.cl,
          ncells = sum(sct_cluster@meta.data$seurat_clusters == cl)
        ),
        10
      )
    }
  )
)
cL_results %>%
  dplyr::group_by(cluster) %>%
  tidyr::nest() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(mm = purrr::map(
    .x = data,
    .f = function(.x) {
      .x %>%
        dplyr::mutate(
          cell = gsub(
            pattern = " \\(.*\\)",
            replacement = "",
            x = type
          )
        ) %>%
        dplyr::mutate(
          tissue = gsub(
            pattern = ".*\\(|\\)",
            replacement = "",
            x = type
          )
        ) %>%
        dplyr::arrange(tissue, -scores)
    }
  ))


sctype_scores <-  cL_results %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 1, wt = scores) %>%
  dplyr::ungroup()
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

sct_cluster@meta.data$sctype <- ""
for(j in unique(sctype_scores$cluster)) {
  cl_type = sctype_scores[sctype_scores$cluster == j, ]
  sct_cluster@meta.data$sctype[sct_cluster@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}



readr::write_rds(
  x = sct_cluster,
  file = "sct_cluster_annotated.rds.gz"
)


p_tsne <- fn_plot_umap(
  .x = sct_cluster,
  .celltype = "sctype",
  .reduction = "tsne"
  )

ggsave(
  filename = "cluster-plot-tsne-sctype.pdf",
  plot = p_tsne,
  device = "pdf",
  width = 13,
  height = 9
)

p_umap <- fn_plot_umap(
  .x = sct_cluster,
  .celltype = "sctype",
  .reduction = "umap"
)

ggsave(
  filename = "cluster-plot-umap-sctype.pdf",
  plot = p_umap,
  device = "pdf",
  width = 10,
  height = 7
)

# write cluster umap ------------------------------------------------------



umap <- as.data.frame(sct_cluster@reductions$umap@cell.embeddings)
colnames(umap) <- c("UMAP_1", "UMAP_2")

cluster <- sct_cluster@meta.data[, c("seurat_clusters", "sctype"), drop=FALSE] %>%
  dplyr::rename(cluster = seurat_clusters)

dplyr::bind_cols(umap, cluster) %>%
  tibble::rownames_to_column(var = "barcode") ->
  cluster_umap

cluster_umap %>%
  readr::write_tsv(
    file = "cluster_umap.tsv",
  )

# Marker genes ------------------------------------------------------------
future::plan(future::sequential)
library(Seurat)
# DefaultAssay(sct_cluster) <- "integrated"
DefaultAssay(sct_cluster) <- "RNA"
sct_cluster <- PrepSCTFindMarkers(object = sct_cluster)

all.markers <- FindAllMarkers(
  object = sct_cluster,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

readr::write_rds(
  x = all.markers,
  file = "sc_sct_cluster_marker_genes.rds.gz"
)

all.markers %>%
  writexl::write_xlsx(
    path = "all.marker.xlsx"
  )

all.markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(n = 3, order_by = avg_log2FC) ->
  all.markers_head

p_marker <- fn_gene_dotplot(
  .sct_cluster = sct_cluster,
  .marker = all.markers,
  .n = 2
)

ggsave(
  filename = "top2-marker-gene-cluster.pdf",
  plot = p_marker,
  device = "pdf",
  width = 13,
  height = 9
)



# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "cellcluster-10x.rda")
