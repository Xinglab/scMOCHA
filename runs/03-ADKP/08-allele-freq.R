#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Aug 24 22:44:15 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
#library(rlang)
library(Seurat)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------
pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
  dplyr::arrange(cancer_types)
  

# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------
fn_load_hetero <- function(.filename) {
  # .filename <- file.path(
  #   sc_dir,
  #   "mgatk_out/final",
  #   "sc.cell_heteroplasmic_df.tsv.gz"
  # )
  
  data.table::fread(input = .filename) |>
    dplyr::rename(barcode = "V1") |>
    tidyr::pivot_longer(
      cols = -barcode,
      names_to = "variant",
      values_to = "af"
    )
}

fn_load_coverage <- function(.filename) {
  
  data.table::fread(
    input = .filename,
    sep = ",",
    col.names = c("pos", "barcode", "depth")
  ) |>
    dplyr::mutate(depth = log2(depth + 1))
}

fn_load_cluster <- function(.filename) {
  data.table::fread(
    input = .filename,
    sep = "\t",
    col.names = c("barcode", "tag", "celltype")
  ) |>
    dplyr::arrange(celltype) |>
    dplyr::mutate(celltype = factor(celltype)) |>
    dplyr::select(-tag)
}

fn_af <- function(.cluster, .hetero) {
  .cluster |>
    dplyr::rename(cluster = celltype) |>
    dplyr::left_join(
      .hetero |> tidyr::pivot_wider(
        names_from  = variant,
        values_from = af
      ),
      by = "barcode"
    )
}

fn_forplot <- function(.af, .coverage) {
  .af |>
    dplyr::select(barcode, cluster, dplyr::contains(">")) |>
    tidyr::pivot_longer(
      cols = -c(barcode, cluster),
      names_to = "variant",
      values_to = "af"
    ) |>
    dplyr::group_by(barcode, cluster) |>
    dplyr::summarise(s_af = sum(af, na.rm = T)) |>
    dplyr::ungroup() |>
    dplyr::arrange(cluster, -s_af) ->
    .rank
  
  .af |>
    dplyr::select(barcode, dplyr::contains(">")) |>
    tidyr::pivot_longer(
      cols = -barcode,
      names_to = "variant",
      values_to = "af"
    ) |>
    dplyr::mutate(
      pos = gsub(pattern = "([[:digit:]]*).*", "\\1", variant) |>
        as.numeric()
    ) |>
    dplyr::left_join(
      .coverage,
      by = c("barcode", "pos")
    ) |>
    tidyr::replace_na(
      replace = list(
        af = 0
      )
    ) |>
    dplyr::mutate(af = ifelse(is.na(depth), NA, af)) |>
    dplyr::arrange(pos) ->
    .forplot
  
  list(
    rank = .rank,
    forplot = .forplot
  )
}

fn_heatmap <- function(.forplot, .cell_variants = NULL, .variant_annotation = NULL) {
  
  .forplot$forplot |>
    dplyr::select(barcode, variant, af)  |>
    tidyr::pivot_wider(
      names_from = "variant",
      values_from = af
    ) |>
    dplyr::slice(
      match(.forplot$rank$barcode, barcode)
    ) |>
    tibble::column_to_rownames(var = "barcode") |>
    as.matrix() |>
    t() ->
    .af_mtx
  
  
  
  tibble::tibble(
    variants = rownames(.af_mtx)
  ) ->
    .for_gcol
  
  .gcol <- if(is.null(.cell_variants)) {
    .for_gcol |>
      dplyr::mutate(
        cell_variants = "black"
      )
  } else {
    .for_gcol |>
      dplyr::mutate(
        cell_variants = ifelse(
          variants %in% .cell_variants,
          "black",
          "red"
        )
      )
  }
  
  
  .forplot$forplot |>
    dplyr::select(barcode, variant, depth) |>
    # dplyr::arrange(pos) |>
    tidyr::pivot_wider(
      names_from = variant,
      values_from = depth
    ) |>
    dplyr::slice(match(.forplot$rank$barcode, barcode)) |>
    tibble::column_to_rownames(var = "barcode") |>
    as.matrix() |>
    t() ->
    .depth_mtx
  
  .forplot$rank |>
    dplyr::select(barcode, cluster) |>
    dplyr::slice(
      match(colnames(.af_mtx), barcode)
    ) |>
    tibble::column_to_rownames(var = "barcode") |>
    dplyr::select(Cluster = cluster) ->
    .af_cluster
  
  
  col_clusters <- levels(.af_cluster$Cluster)
  col_colors <- pcc$color[1:length(levels(.af_cluster$Cluster))]
  
  names(col_colors) <- col_clusters
  
  chm_top <- ComplexHeatmap::HeatmapAnnotation(
    df = .af_cluster,
    gap = unit(c(2,2), "mm"),
    col = list(Cluster = col_colors),
    which = "column"
  )
  
  ch_af <- if(!is.null(.variant_annotation)) {
    
    hma_right <- ComplexHeatmap::rowAnnotation(
      df = .variant_annotation |>
        dplyr::select(-Conservation)
    )
    hma_left <- ComplexHeatmap::rowAnnotation(
      df = .variant_annotation |>
        dplyr::select(Conservation)
    )
    ComplexHeatmap::Heatmap(
      matrix = .af_mtx,
      col = circlize::colorRamp2(
        breaks = c(0, 0.98, 1),
        colors = c("white", "red", "#440154FF"),
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
      row_names_gp = gpar(
        # fontsize = 20,
        col = .gcol$cell_variants
      ),
      # column
      cluster_columns = FALSE,
      cluster_column_slices = T,
      # clustering_distance_columns = "pearson",
      # clustering_method_columns = "ward.D",
      show_column_names = FALSE,
      row_names_side = "left",
      
      top_annotation = chm_top,
      left_annotation = hma_left,
      right_annotation = hma_right
    )
    
  } else {
    ComplexHeatmap::Heatmap(
      matrix = .af_mtx,
      col = circlize::colorRamp2(
        breaks = c(0, 0.98, 1),
        colors = c("white", "red", "#440154FF"),
        space = "RGB"
      ),
      name = "Allele Freq",
      na_col = "grey",
      color_space = "LAB",
      # rect_gp = gpar(col = NA),
      border = NA,
      cell_fun = NULL,
      layer_fun = NULL,
      jitter = FALSE,
      # row
      cluster_rows = F,
      cluster_row_slices = T,
      clustering_distance_rows = "pearson",
      clustering_method_rows = "ward.D",
      # row_names_gp = gpar(
      #   col = .gcol$cell_variants
      # ),
      # column
      cluster_columns = FALSE,
      cluster_column_slices = T,
      # clustering_distance_columns = "pearson",
      # clustering_method_columns = "ward.D",
      show_column_names = FALSE,
      row_names_side = "left",
      
      top_annotation = chm_top,
    )
    
  }
  
  
  ComplexHeatmap::Heatmap(
    matrix = .depth_mtx,
    col = circlize::colorRamp2(
      breaks = c(0, quantile(.depth_mtx, na.rm = T, probs = 0.75)),
      colors = c("white", "red"),
      # colors =  c("#440154FF", "#FDE725FF"),
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
    row_names_gp = gpar(
      # fontsize = 20,
      col = .gcol$cell_variants
    ),
    # column
    cluster_columns = FALSE,
    cluster_column_slices = T,
    # clustering_distance_columns = "pearson",
    # clustering_method_columns = "ward.D",
    show_column_names = FALSE,
    row_names_side = "left",
    
    top_annotation = chm_top
  ) ->
    ch_depth
  
  list(
    ch_af = ch_af,
    ch_depth = ch_depth
  )
}

# load data ---------------------------------------------------------------

metadata_anno_depth <- readr::read_rds(
  file = "/home/liuc9/github/scMOCHA/03-ADKP/output/rda/metadata_anno_depth.rds"
)

# body --------------------------------------------------------------------
colnames(metadata_anno_depth)


metadata_anno_depth |>
  dplyr::filter(!is.na(linkfile)) |> 
  dplyr::mutate(dia = `Diagnosis (neurology)`) |> 
  dplyr::mutate(
    source_name = dia
  ) |> 
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
  dplyr::mutate(
    Sex = factor(Sex),
    dia = factor(dia)
  ) |> 
  tidyr::unnest(cols = n_na) |> 
  dplyr::filter(
    dep_med > 1000
  ) |> 
  dplyr::mutate(
    variant = purrr::map2(
      .x = anno,
      .y = tardir,
      .f = function(.x, .y) {
        if(is.na(.y)) {return(NULL)}
        .x |>
          dplyr::mutate(
            variant = glue::glue("{Position}{Ref}>{Alt}")
          ) |>
          dplyr::select(variant)
      }
    )
  ) |> 
  dplyr::mutate(srrid = `Sample ID`) ->
  metadata_anno_depth_variant

metadata_anno_depth_variant |> 
  dplyr::mutate(
    color = dplyr::case_match(
      source_name,
      "AD" ~ggsci::pal_jama()(4)[[2]],
      "MCI" ~ ggsci::pal_jama()(4)[[1]],
      )
    ) |>
  dplyr::select(srrid, source_name, variant, color) |>
  dplyr::filter(!purrr::map_lgl(.x = variant, .f = is.null)) ->
  for_variant

fn_upset_plot <- function(.x) {
  # .x <- "nCoV_PBMC(severe)"
  library(ggupset)
  for_variant |>
    dplyr::filter(source_name == .x) ->
    d
  
  d |>
    tidyr::unnest(cols = variant) |>
    dplyr::select(-source_name) |>
    dplyr::group_by(variant) |>
    tidyr::nest() |>
    dplyr::ungroup() |>
    dplyr::mutate(
      srrid = purrr::map(
        .x = data,
        .f = function(.x) {
          .x |> dplyr::pull(srrid)
        }
      )
    ) |>
    dplyr::select(-data) ->
    dd
  
  
  dd |> 
    ggplot(aes(x = srrid)) +
    geom_bar(width = 0.6, fill = d$color[1]) +
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
      combmatrix.panel.point.color.fill = d$color[1],
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
      title = .x
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(size = 0.5, color = "black"),
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
    .p_up
  
  ggsave(
    plot = .p_up,
    filename = "upset-{.x}.pdf" |> glue::glue(),
    path = "/home/liuc9/github/scMOCHA/03-ADKP/output",
    width = 9,
    height = 6,
    device = "pdf"
  )
  
  dd |> 
    dplyr::mutate(n = purrr::map_int(
      .x = srrid, 
      .f = length
    )) |> 
    dplyr::mutate(
      sharing = purrr::map_chr(
        .x = srrid,
        .f = paste0,
        collapse = ","
      )
    ) |> 
    dplyr::arrange(n) |> 
    dplyr::select(-srrid) -> 
    .v
  
  list(
    v = .v,
    p_up = .p_up
  )
  
}

metadata_anno_depth_variant$source_name |>
  unique() |>
  purrr::map(
    .f = fn_upset_plot
  ) ->
  p_ups

(p_ups[[1]]$p_up / p_ups[[2]]$p_up) +
  plot_annotation(tag_levels = "A") ->
  p_ups_together;p_ups_together

ggsave(
  plot = p_ups_together,
  filename = "upset-all.pdf" |> glue::glue(),
  path = "/home/liuc9/github/scMOCHA/03-ADKP/output",
  width = 12,
  height = 10,
  device = "pdf"
)


# save to xlsx ------------------------------------------------------------


names(p_ups) <- unique(metadata_anno_depth_variant$source_name)

p_ups |> 
  purrr::map("v") |> 
  writexl::write_xlsx(
    path = "/home/liuc9/github/scMOCHA/03-ADKP/output/upset-variants.xlsx"
  )



for_variant |>
  dplyr::select(srrid, variant) |>
  dplyr::mutate(
    variant = purrr::map(
      .x = variant,
      .f = function(.x) {
        .x |> dplyr::pull(variant)
      }
    )
  ) |> 
  dplyr::filter(grepl(
    pattern = "AD",
    x = srrid,
  )) |> 
  tibble::deframe() |> 
  purrr::reduce(.f = intersect) -> 
  common_variants_AD

for_variant |>
  dplyr::select(srrid, variant) |>
  dplyr::mutate(
    variant = purrr::map(
      .x = variant,
      .f = function(.x) {
        .x |> dplyr::pull(variant)
      }
    )
  ) |> 
  dplyr::filter(grepl(
    pattern = "MCI",
    x = srrid,
  )) |> 
  tibble::deframe() |> 
  purrr::reduce(.f = intersect) -> 
  common_variants_MCI


p_text_mci <- ggplot() +
  annotate(
    "text", x = 1, y = 1,
    size = 6,
    color = "black", #ggsci::pal_aaas()(1),
    label = stringr::str_wrap(
      string = common_variants_MCI |> 
        paste0(collapse = ", "),
      width = 30
    )
  ) +
  theme_void()

p_text_ad <- ggplot() +
  annotate(
    "text", x = 1, y = 1,
    size = 6,
    color = "black", #ggsci::pal_aaas()(1),
    label = stringr::str_wrap(
      string = common_variants_AD |> 
        paste0(collapse = ", "),
      width = 30
    )
  ) +
  theme_void()

p_text_ad / p_text_mci

# All variants ------------------------------------------------------------



for_variant |>
  dplyr::select(srrid, variant) |>
  dplyr::mutate(
    variant = purrr::map(
      .x = variant,
      .f = function(.x) {
        .x |> dplyr::pull(variant)
      }
    )
  ) |> 
  tibble::deframe() |> 
  purrr::reduce(.f = union) ->
  all_variants


tibble::tibble(
  variant = all_variants,
  pos = gsub(
    pattern = ">|[A-Z]",
    replacement = "",
    x = variant
  )
) |> 
  dplyr::mutate(pos = as.integer(pos)) ->
  all_variants_pos


future::plan(future::multisession, workers = 10)
metadata_anno_depth_variant |> 
  dplyr::select(srrid, tardir) |> 
  dplyr::mutate(
    a = furrr::future_map2(
      .x = tardir,
      .y = srrid,
      .f = \(.tardir, .srrid) {
        
        cluster_umap <- fn_load_cluster(
          .filename = file.path(
            .tardir,
            "barcode_cluster.tsv"
          )
        ) |> 
          dplyr::mutate(barcode = glue::glue("{.srrid}#{barcode}"))
        
        cell_hetero_raw <- fn_load_hetero(
          .filename = file.path(
            .tardir,
            "cell.cell_heteroplasmic_df_raw.tsv.gz"
          )
        ) |> 
          dplyr::filter(variant %in% all_variants_pos$variant) |> 
          dplyr::mutate(barcode = glue::glue("{.srrid}#{barcode}")) |> 
          dplyr::left_join(
            cluster_umap,
            by = "barcode"
          ) |> 
          dplyr::rename(cluster = celltype)
        
        cell_coverage <- fn_load_coverage(
          .filename = file.path(
            .tardir,
            "cell.coverage.txt.gz"
          )
        ) |> 
          dplyr::filter(pos %in% all_variants_pos$pos) |> 
          dplyr::mutate(barcode = glue::glue("{.srrid}#{barcode}"))
        
        tibble::tibble(
          cell_hetero_raw = list(cell_hetero_raw),
          cell_coverage = list(cell_coverage)
        )
        
      }
    )
  ) ->
  metadata_anno_depth_variant_hetero_coverage
future::plan(future::sequential)


af <- metadata_anno_depth_variant_hetero_coverage |> 
  dplyr::select(srrid, a) |> 
  tidyr::unnest(cols = a) |> 
  dplyr::select(cell_hetero_raw) |> 
  tidyr::unnest(cols = cell_hetero_raw) |> 
  tidyr::pivot_wider(
    names_from = variant,
    values_from = af
  ) |> 
  dplyr::filter(!is.na(cluster))

covera <- metadata_anno_depth_variant_hetero_coverage |> 
  dplyr::select(srrid, a) |> 
  tidyr::unnest(cols = a) |>  
  dplyr::select(cell_coverage) |> 
  tidyr::unnest(cols = cell_coverage)

cell_raw_cluster_forplot <- fn_forplot(
  .af = af,
  .coverage = covera
)

pp <- fn_heatmap(
  .forplot = cell_raw_cluster_forplot
)



# Cell cluster ------------------------------------------------------------
sctu_tsne <- readr::read_rds("/home/liuc9/github/scMOCHA/03-ADKP/forrefs/azimuth_syn21438358/sctu_tsne.rds")

sctu_tsne$seurat_clusters |> 
  table() ->
  cl

ccl <- factor(glue::glue("{names(cl)} (n={cl})"), glue::glue("{names(cl)} (n={cl})"))

sctu_tsne$seurat_clusters_n <- ccl[sctu_tsne$seurat_clusters]



DimPlot(
  object = sctu_tsne,
  reduction = "tsne",
  cols = paletteer::paletteer_d(
    palette = "ggsci::springfield_simpsons",
    direction = -1
  ),
  group.by = "seurat_clusters_n"
)

sctu_tsne@reductions$tsne@cell.embeddings

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------

# save.image(file = "/scr1/users/liuc9/mitochondrial/realdata/03-ADKP/output/rda/08-allele-freq.rda")

load(file = "/scr1/users/liuc9/mitochondrial/realdata/03-ADKP/output/rda/08-allele-freq.rda")
