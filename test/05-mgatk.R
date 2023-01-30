# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Jan 26 00:13:25 2023
# @DESCRIPTION: 05-mgatk.R

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
hetero_file <- "/home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/MTbam/mgatk/final/mgatk.cell_heteroplasmic_df.tsv.gz"

hetero <- vroom::vroom(file = hetero_file) %>% 
  dplyr::rename(barcode = `...1`) %>% 
  tidyr::pivot_longer(
    cols = -barcode,
    names_to = "variant",
    values_to = "af"
  ) 


# body --------------------------------------------------------------------

hetero %>% 
  dplyr::group_by(variant) %>% 
  dplyr::summarise(s_af = sum(af)) %>% 
  dplyr::arrange(s_af) ->
  hetero_variant_rank

hetero %>% 
  dplyr::group_by(barcode) %>% 
  dplyr::summarise(s_af = sum(af)) %>% 
  dplyr::arrange(s_af) ->
  hetero_barcode_rank

hetero %>% 
  dplyr::mutate(variant = factor(variant, levels = hetero_variant_rank$variant)) %>% 
  dplyr::mutate(barcode = factor(barcode, levels = hetero_barcode_rank$barcode)) %>% 
  dplyr::rename(`Allele Freq` = af) %>% 
  ggplot(aes(
    x = barcode,
    y = variant
  )) +
  geom_tile(aes(fill = `Allele Freq`)) +
  # scale_fill_gradient2(
  #   low = "white",
  #   mid = CPCOLS[3],
  #   high = "#67000e"
  # ) +
  # scale_fill_brewer(name = "Alelle Freq") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    x = "Cell",
    y = "MT variant"
  ) ->
  tileplot;tileplot

ggsave(
  filename = "mgatk_cell_genotype_tileplot.pdf",
  plot = tileplot,
  device = "pdf",
  path = "data/PBMC_10k_v3_10x/result/04-allele-freq/",
  width = 12,
  height = 9
)

depthtable <- readr::read_tsv(file = "/home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/MTbam/mgatk/final/mgatk.depthTable.txt", col_names = c("barcode", "depth"))

depthtable %>% 
  ggplot(aes(x = depth)) +
  geom_density()

depthtable$depth %>% quantile()


coverage <- vroom::vroom(file = "/scr1/users/liuc9/mitochondrial/PBMC_10k_v3_10x/MTbam/mgatk/final/mgatk.coverage.txt.gz", delim = ",", col_names = c("pos", "barcode", "depth"))

coverage %>% 
  dplyr::mutate(depth = log2(depth + 1)) ->
  coverage_log2

coverage_log2 %>% 
  tidyr::pivot_wider(
    names_from = pos,
    values_from = depth
  ) ->
  coverage_wider

# seurat ------------------------------------------------------------------



sct_cluster <- readr::read_rds("data/PBMC_10k_v3_10x/rda/pbmc_sct_cluster_annotated.rds.gz")


umap <- sct_cluster@reductions$umap@cell.embeddings
colnames(umap) <- c("UMAP_1", "UMAP_2")

cluster <- sct_cluster@meta.data[, c("seurat_clusters", "sctype")] %>% 
  dplyr::rename(cluster = seurat_clusters, celltype =  sctype)

hetero %>% 
  tidyr::pivot_wider(names_from  = variant, values_from = af) ->
  hetero_w

dplyr::bind_cols(umap, cluster) %>% 
  tibble::rownames_to_column(var = "barcode") %>% 
  dplyr::left_join(hetero_w, by = "barcode") ->
  # dplyr::mutate(
  #   dplyr::across(
  #   dplyr::everything(), 
  #   ~tidyr::replace_na(.x, 0)
  #   )) ->
  cell_cluster_af

cell_cluster_af %>% 
  dplyr::left_join(depthtable, by = "barcode") ->
  cell_cluster_af_depth

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
myPalette <- colorRampPalette(c("#969696", "#fa0202"))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 1))

cell_cluster_af %>% 
  # dplyr::filter(variant == "2617A>G") %>% nrow()
  ggplot() +
  geom_point(
    aes(
      x = UMAP_1,
      y = UMAP_2,
      colour = `9966G>A`,
      shape = NULL,
      alpha = NULL
    ),
    size = 0.7
  ) +
  sc +
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
  ) ->
  p_variant;p_variant


ggsave(
  filename = "mgatk_cell_genotype_variant.pdf",
  plot = p_variant,
  device = "pdf",
  path = "data/PBMC_10k_v3_10x/result/04-allele-freq/",
  width = 12,
  height = 9
)


cell_cluster_af %>% 
  # dplyr::filter(variant == "2617A>G") %>% nrow()
  ggplot() +
  geom_point(
    aes(
      x = UMAP_1,
      y = UMAP_2,
      colour = `263A>G`,
      shape = NULL,
      alpha = NULL
    ),
    size = 0.7
  ) +
  sc +
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
  ) ->
  p_variant_ag;p_variant_ag


ggsave(
  filename = "mgatk_cell_genotype_variant_ag.pdf",
  plot = p_variant_ag,
  device = "pdf",
  path = "data/PBMC_10k_v3_10x/result/04-allele-freq/",
  width = 12,
  height = 9
)


cell_cluster_af_depth %>% 
  ggplot(aes(
    x = depth,
    y = `9966G>A`,
  )) +
  geom_point() +
  scale_x_continuous(
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    expand = c(0.01, 0)
  ) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(
      colour = "black",
      size = 0.5,
      # arrow = grid::arrow(
      #   angle = 1,
      #   length = unit(1, "npc"),
      #   type = "closed"
      # )
    ),
    # axis.ticks = element_blank(),
    axis.text = element_text(
      size = 12,
      color = "black",
      face = "bold"
    ),
    axis.title = element_text(
      size = 14, 
      face = "bold", 
      # hjust = 0.01
    ),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(
      face = "bold",
      color = "black",
      size = 12
    )
  ) +
  labs(
    x = "Average depth"
  )


cell_cluster_af_depth %>% 
  dplyr::select(barcode, depth, `9966G>A`) %>% 
  dplyr::inner_join(
    coverage_wider %>% 
      dplyr::select(barcode, `9966`),
    by = "barcode"
  ) %>% 
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(), 
      ~tidyr::replace_na(.x, 0)
    )) ->
  cell_cluster_af_depth_select

cor.test(cell_cluster_af_depth_select$`9966G>A`, cell_cluster_af_depth_select$`9966`) -> pearson

cell_cluster_af_depth_select %>% 
  ggplot(aes(
    x = `9966`,
    y = `9966G>A`,
  )) +
  geom_point(position = position_jitter()) +
  scale_x_continuous(
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    expand = c(0.01, 0)
  ) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(
      colour = "black",
      size = 0.5,
      # arrow = grid::arrow(
      #   angle = 1,
      #   length = unit(1, "npc"),
      #   type = "closed"
      # )
    ),
    # axis.ticks = element_blank(),
    axis.text = element_text(
      size = 12,
      color = "black",
      face = "bold"
    ),
    axis.title = element_text(
      size = 14, 
      face = "bold", 
      # hjust = 0.01
    ),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(
      face = "bold",
      color = "black",
      size = 12
    )
  ) +
  annotate(
    geom = "text",
    size = 11,
    x = 40,
    y = 0.8,
    label = glue::glue(
      "$Pearson \\ r=<<round(pearson$estimate,3)>>$",
      .open = "<<", 
      .close = ">>"
    ) %>% 
      latex2exp::TeX(),
  ) +
  labs(
    x = "chrMT:9966 read depth",
    y = "9966G>A Allele Freq"
  ) ->
  p_variant_cor;p_variant_cor


ggsave(
  filename = "mgatk_cell_genotype_variant_cor.pdf",
  plot = p_variant_cor,
  device = "pdf",
  path = "data/PBMC_10k_v3_10x/result/04-allele-freq/",
  width = 7,
  height = 6
)




cell_cluster_af_depth %>% 
  dplyr::select(barcode, depth, `263A>G`) %>% 
  dplyr::inner_join(
    coverage_wider %>% 
      dplyr::select(barcode, `263`),
    by = "barcode"
  ) %>% 
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(), 
      ~tidyr::replace_na(.x, 0)
    )) ->
  cell_cluster_af_depth_select

cor.test(cell_cluster_af_depth_select$`263A>G`, cell_cluster_af_depth_select$`263`) -> pearson

cell_cluster_af_depth_select %>% 
  ggplot(aes(
    x = `263`,
    y = `263A>G`,
  )) +
  # geom_point() +
  geom_point(position = position_jitter(
    width = 0.5,
    height = 0.002
  )) +
  scale_x_continuous(
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    expand = c(0.01, 0)
  ) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(
      colour = "black",
      size = 0.5,
      # arrow = grid::arrow(
      #   angle = 1,
      #   length = unit(1, "npc"),
      #   type = "closed"
      # )
    ),
    # axis.ticks = element_blank(),
    axis.text = element_text(
      size = 12,
      color = "black",
      face = "bold"
    ),
    axis.title = element_text(
      size = 14, 
      face = "bold", 
      # hjust = 0.01
    ),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(
      face = "bold",
      color = "black",
      size = 12
    )
  ) +
  annotate(
    geom = "text",
    size = 8,
    x = 2,
    y = 0.8,
    label = glue::glue(
      "$Pearson \\ r=<<round(pearson$estimate,3)>>$",
      .open = "<<", 
      .close = ">>"
    ) %>% 
      latex2exp::TeX(),
  ) +
  labs(
    x = "chrMT:263 read depth",
    y = "263A>G Allele Freq"
  ) ->
  p_variant_cor;p_variant_cor


ggsave(
  filename = "mgatk_cell_genotype_variant_cor_ag.pdf",
  plot = p_variant_cor,
  device = "pdf",
  path = "data/PBMC_10k_v3_10x/result/04-allele-freq/",
  width = 7,
  height = 6
)

# 
# cell_cluster_af %>% 
#   dplyr::select(barcode, dplyr::contains(">")) %>% 
#   tidyr::pivot_longer(
#     cols = -barcode,
#     names_to = "variant",
#     values_to = "af"
#   )




# Heatmap -----------------------------------------------------------------



cell_cluster_af %>% 
  dplyr::select(barcode, dplyr::contains(">")) %>% 
  tidyr::pivot_longer(
    cols = -barcode,
    names_to = "variant",
    values_to = "af"
  ) %>% 
  dplyr::group_by(barcode) %>% 
  dplyr::summarise(s_af = sum(af)) %>% 
  dplyr::left_join(
    cell_cluster_af %>%
      dplyr::select(barcode, cluster) ,
    by = "barcode"
  ) %>% 
  dplyr::arrange(cluster, -s_af) ->
  cell_cluster_af_col_rank

cell_cluster_af %>% 
  dplyr::select(barcode, dplyr::contains(">")) %>% 
  tidyr::pivot_longer(
    cols = -barcode,
    names_to = "variant",
    values_to = "af"
  ) %>% 
  dplyr::mutate(pos = gsub(pattern = "([[:digit:]]*).*", "\\1", variant) %>% as.numeric()) ->
  cell_cluster_af_pos


coverage_log2 %>% 
  dplyr::filter(barcode %in% cell_cluster_af_pos$barcode) %>% 
  dplyr::filter(pos %in% cell_cluster_af_pos$pos) ->
  coverage_log2_pos

cell_cluster_af_pos %>% 
  dplyr::left_join(
    coverage_log2_pos,
    by = c("barcode", "pos")
  ) %>% 
  dplyr::mutate(af = ifelse(is.na(depth), NA, af)) %>% 
  dplyr::arrange(pos) %>% 
  dplyr::select(barcode, variant, af) %>% 
  tidyr::pivot_wider(
    names_from = "variant",
    values_from = af
  ) %>% 
  dplyr::slice(match(cell_cluster_af_col_rank$barcode, barcode)) %>% 
  tibble::column_to_rownames(var = "barcode") %>% 
  # dplyr::filter_all(
  #   dplyr::any_vars(.!=0)
  # ) %>% 
  as.matrix() %>% 
  t() ->
  cell_cluster_af_mtx 

library(ComplexHeatmap)
pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv")


cell_cluster_af %>% 
  dplyr::select(barcode, cluster, celltype ) %>% 
  dplyr::slice(match(colnames(cell_cluster_af_mtx), barcode)) %>% 
  tibble::column_to_rownames(var = "barcode") %>% 
  # dplyr::select(celltype) %>% 
  # dplyr::rename("Cell type" = celltype) ->
  dplyr::select(Cluster = cluster) ->
  cell_cluster_af_cluster

col_clusters <- levels(cell_cluster_af_cluster$Cluster) %>% as.numeric()
col_colors <- pcc$color[1:length(levels(cell_cluster_af_cluster$Cluster))]

names(col_colors) <- col_clusters

chm_top <- ComplexHeatmap::HeatmapAnnotation(
  df = cell_cluster_af_cluster,
  gap = unit(c(2,2), "mm"),
  col = list(Cluster = col_colors),
  which = "column"
)


ComplexHeatmap::Heatmap(
  matrix = cell_cluster_af_mtx,
  col = circlize::colorRamp2(
    breaks = c(0, 1), 
    colors = c("white", "red"), 
    space = "RGB"
  ),
  name = "Allele Freq",
  na_col = "grey",
  color_space = "LAB",
  rect_gp = gpar(col = NA),
  border = NA,
  cell_fun = NULL,
  layer_fun = NULL,
  jitter = FALSE,
  # row
  cluster_rows = F,
  cluster_row_slices = T,
  clustering_distance_rows = "pearson",
  clustering_method_rows = "ward.D",
  # column
  cluster_columns = FALSE,
  cluster_column_slices = T,
  # clustering_distance_columns = "pearson",
  # clustering_method_columns = "ward.D",
  show_column_names = FALSE,
  
  top_annotation = chm_top
) ->
  chm;chm


{
  pdf(
    file = "data/PBMC_10k_v3_10x/result/04-allele-freq/mgatk_cell_genotype_heatmap.pdf",
    width = 12, 
    height = 4
  )
  ComplexHeatmap::draw(object = chm)
  
  dev.off()
}

cell_cluster_af_pos %>% 
  dplyr::left_join(
    coverage_log2_pos,
    by = c("barcode", "pos")
  ) %>% 
  dplyr::select(barcode, pos, depth) %>% 
  dplyr::arrange(pos) %>% 
  tidyr::pivot_wider(
    names_from = pos,
    values_from = depth
  ) %>% 
  dplyr::slice(match(cell_cluster_af_col_rank$barcode, barcode)) %>% 
  tibble::column_to_rownames(var = "barcode") %>% 
  # dplyr::filter_all(
  #   dplyr::any_vars(.!=0)
  # ) %>% 
  as.matrix() %>% 
  t() ->
  depth_mtx



ComplexHeatmap::Heatmap(
  matrix = depth_mtx,
  col = circlize::colorRamp2(
    breaks = c(1, 4), 
    colors = c("white", "red"), 
    space = "RGB"
  ),
  name = "log2(Depth+1)",
  na_col = "grey",
  color_space = "LAB",
  rect_gp = gpar(col = NA),
  border = NA,
  cell_fun = NULL,
  layer_fun = NULL,
  jitter = FALSE,
  # row
  cluster_rows = F,
  cluster_row_slices = T,
  clustering_distance_rows = "pearson",
  clustering_method_rows = "ward.D",
  # column
  cluster_columns = FALSE,
  cluster_column_slices = T,
  # clustering_distance_columns = "pearson",
  # clustering_method_columns = "ward.D",
  show_column_names = FALSE,
  
  top_annotation = chm_top
) ->
  chm_depth;chm_depth

{
  pdf(
    file = "data/PBMC_10k_v3_10x/result/04-allele-freq/mgatk_cell_depth_heatmap.pdf",
    width = 12, 
    height = 4
  )
  ComplexHeatmap::draw(object = chm_depth)
  
  dev.off()
}


coverage_log2


coverage_log2 %>% 
  # dplyr::arrange(pos) %>% 
  tidyr::pivot_wider(
    names_from = pos,
    values_from = depth
  ) %>% 
  dplyr::slice(match(cell_cluster_af_col_rank$barcode, barcode)) %>% 
  tibble::column_to_rownames(var = "barcode") %>% 
  # dplyr::filter_all(
  #   dplyr::any_vars(.!=0)
  # ) %>% 
  as.matrix() %>% 
  t() ->
  depth_all_mtx



ComplexHeatmap::Heatmap(
  matrix = depth_all_mtx,
  col = circlize::colorRamp2(
    breaks = c(1, 10), 
    colors = c("white", "red"), 
    space = "RGB"
  ),
  name = "log2(Depth+1)",
  na_col = "grey",
  color_space = "LAB",
  rect_gp = gpar(col = NA),
  border = NA,
  cell_fun = NULL,
  layer_fun = NULL,
  jitter = FALSE,
  # row
  cluster_rows = F,
  cluster_row_slices = T,
  clustering_distance_rows = "pearson",
  clustering_method_rows = "ward.D",
  # column
  cluster_columns = FALSE,
  cluster_column_slices = T,
  # clustering_distance_columns = "pearson",
  # clustering_method_columns = "ward.D",
  show_column_names = FALSE,
  show_row_names = FALSE,
  
  top_annotation = chm_top
) ->
  chm_depth_all;chm_depth_all

{
  pdf(
    file = "data/PBMC_10k_v3_10x/result/04-allele-freq/mgatk_cell_depth_heatmap.pdf",
    width = 12, 
    height = 30
  )
  ComplexHeatmap::draw(object = chm_depth)
  
  dev.off()
}


# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "data/PBMC_10k_v3_10x/rda/05-mgatk.rda")

load(file = "data/PBMC_10k_v3_10x/rda/05-mgatk.rda")
