# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Fri Mar 24 00:44:44 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |> 
  dplyr::arrange(cancer_types)

# src ---------------------------------------------------------------------


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
    sep = "\t"
  ) |> 
    dplyr::mutate(cluster = factor(cluster)) |> 
    dplyr::mutate(sctype = glue::glue("{cluster}, {sctype}")) |> 
    dplyr::mutate(sctype = factor(sctype))
}

fn_af <- function(.cluster, .hetero) {
  .cluster |>
    dplyr::select(-cluster) |> 
    dplyr::rename(cluster = sctype) |> 
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
# load data ---------------------------------------------------------------
sc_dir <- "/home/liuc9/scratch/mitochondrial/testdata/1_old_donor_pbmc/flu2/Flu2/outs"
flu2_cluster_umap <- fn_load_cluster(
  .filename = file.path(
    sc_dir,
    "cluster_umap.tsv"
  )
)
flu2_coverage <- fn_load_coverage(
  .filename = file.path(
    sc_dir,
    "mgatk_out/final",
    "sc.coverage.txt.gz"
  ) 
)
flu2_hetero_raw <- fn_load_hetero(
  .filename = file.path(
    sc_dir,
    "mgatk_out/final",
    "sc.cell_heteroplasmic_df_raw.tsv.gz"
  )
)
flu2_hetero_cluster <- fn_load_hetero(
  .filename = file.path(
    sc_dir,
    "mgatk_cluster/final",
    "mgatk_cluster.cell_heteroplasmic_df.tsv.gz"
  )
) |> 
  dplyr::mutate(
    barcode = gsub(
      pattern = "cluster|-1",
      replacement = "",
      x = barcode
    )
  ) |> 
  dplyr::mutate(barcode = as.integer(barcode) -1) |> 
  dplyr::mutate(cluster = barcode) |> 
  dplyr::mutate(cluster = factor(cluster)) |> 
  dplyr::left_join(
    flu2_cluster_umap |> 
      dplyr::select(cluster, sctype) |> 
      dplyr::distinct(),
    by = "cluster"
  ) |> 
  dplyr::select(-cluster) |> 
  dplyr::rename(cluster = sctype)
flu2_cell_raw_cluster_af <- flu2_cluster_umap |> 
  dplyr::left_join(flu2_hetero_raw, by = "barcode") |> 
  dplyr::select(-cluster) |> 
  dplyr::rename(cluster = sctype) |> 
  dplyr::filter(variant %in% flu2_hetero_cluster$variant) |> 
  tidyr::pivot_wider(
    names_from = variant,
    values_from = af
  )
flu2_cell_raw_cluster_forplot <- fn_forplot(
  .af = flu2_cell_raw_cluster_af, 
  .coverage = flu2_coverage
)





sc_dir <- "/home/liuc9/scratch/mitochondrial/testdata/1_old_donor_pbmc/flu5/Flu5/outs"
flu5_cluster_umap <- fn_load_cluster(
  .filename = file.path(
    sc_dir,
    "cluster_umap.tsv"
  )
)
flu5_coverage <- fn_load_coverage(
  .filename = file.path(
    sc_dir,
    "mgatk_out/final",
    "sc.coverage.txt.gz"
  ) 
)
flu5_hetero_raw <- fn_load_hetero(
  .filename = file.path(
    sc_dir,
    "mgatk_out/final",
    "sc.cell_heteroplasmic_df_raw.tsv.gz"
  )
)
flu5_hetero_cluster <- fn_load_hetero(
  .filename = file.path(
    sc_dir,
    "mgatk_cluster/final",
    "mgatk_cluster.cell_heteroplasmic_df.tsv.gz"
  )
) |> 
  dplyr::mutate(
    barcode = gsub(
      pattern = "cluster|-1",
      replacement = "",
      x = barcode
    )
  ) |> 
  dplyr::mutate(barcode = as.integer(barcode) -1) |> 
  dplyr::mutate(cluster = barcode) |> 
  dplyr::mutate(cluster = factor(cluster)) |> 
  dplyr::left_join(
    flu5_cluster_umap |> 
      dplyr::select(cluster, sctype) |> 
      dplyr::distinct(),
    by = "cluster"
  ) |> 
  dplyr::select(-cluster) |> 
  dplyr::rename(cluster = sctype)
flu5_cell_raw_cluster_af <- flu5_cluster_umap |> 
  dplyr::left_join(flu5_hetero_raw, by = "barcode") |> 
  dplyr::select(-cluster) |> 
  dplyr::rename(cluster = sctype) |> 
  dplyr::filter(variant %in% flu5_hetero_cluster$variant) |> 
  tidyr::pivot_wider(
    names_from = variant,
    values_from = af
  )
flu5_cell_raw_cluster_forplot <- fn_forplot(
  .af = flu5_cell_raw_cluster_af, 
  .coverage = flu5_coverage
)

# body --------------------------------------------------------------------


flu2_cell_raw_cluster_forplot$forplot |> 
  dplyr::select(barcode, variant, af, depth) |> 
  dplyr::left_join(
    flu2_cell_raw_cluster_forplot$rank |> 
      dplyr::select(barcode, cluster),
    by = "barcode"
  ) |> 
  tidyr::separate_wider_delim(
    cols = cluster,
    delim = ", ",
    names = c("cluster", "sctype")
  ) ->
  flu2_variant
  


flu5_cell_raw_cluster_forplot$forplot |> 
  dplyr::select(barcode, variant, af, depth) |> 
  dplyr::left_join(
    flu5_cell_raw_cluster_forplot$rank |> 
      dplyr::select(barcode, cluster),
    by = "barcode"
  ) |> 
  tidyr::separate_wider_delim(
    cols = cluster,
    delim = ", ",
    names = c("cluster", "sctype")
  ) ->
  flu5_variant


# All variants ------------------------------------------------------------


ggvenn::ggvenn(
  data = list(
    Flu2 = unique(flu2_variant$variant),
    Flu5 = unique(flu5_variant$variant)
  ),
  show_percentage = FALSE,
  fill_color = ggsci::pal_aaas()(3),
  fill_alpha = 0.9,
  stroke_size = 0.5,
  set_name_size = 10,
  text_size = 8
)



ggplot() +  
  annotate(
    "text", x = 1, y = 1,
    size = 6,
    label = stringr::str_wrap(
      string = intersect(
        unique(flu2_variant$variant),
        unique(flu5_variant$variant)
      ) |> paste0(collapse = ", "),
      width = 30
    )
    ) + 
  theme_void()

intersect(unique(flu2_variant$variant), unique(flu5_variant$variant)) |> 
  paste0(collapse = ",")

ggplot() +  
  annotate(
    "text", x = 1, y = 1,
    size = 6,
    color = ggsci::pal_aaas()(1),
    label = stringr::str_wrap(
      string = setdiff(
        unique(flu2_variant$variant),
        unique(flu5_variant$variant)
      ) |> paste0(collapse = ", "),
      width = 30
    )
  ) + 
  theme_void()

setdiff(
  unique(flu2_variant$variant),
  unique(flu5_variant$variant)
) |> paste0(collapse = ", ")

ggplot() +  
  annotate(
    "text", x = 1, y = 1,
    size = 6,
    color = ggsci::pal_aaas()(2)[2],
    label = stringr::str_wrap(
      string = setdiff(
        unique(flu5_variant$variant),
        unique(flu2_variant$variant)
      ) |> paste0(collapse = ", "),
      width = 30
    )
  ) + 
  theme_void()

setdiff(
  unique(flu5_variant$variant),
  unique(flu2_variant$variant)
) |> paste0(collapse = ", ")

# cluster -----------------------------------------------------------------

inter_celltype <- intersect(
  unique(flu2_variant$sctype),
  unique(flu5_variant$sctype)
)

flu2_variant |> 
  dplyr::filter(sctype %in% inter_celltype) |> 
  dplyr::select(variant, af, sctype) |> 
  dplyr::group_by(variant) |> 
  tidyr::nest() |> 
  dplyr::ungroup() |> 
  dplyr::mutate(
    .ratio = purrr::map_dbl(
      .x = data,
      .f = function(.x) {
        .ncell <- nrow(.x)
        .x |> 
          dplyr::filter(!is.na(af)) |> 
          dplyr::filter(af != 0) |> 
          nrow() ->
          .n_detected
        
        .ratio <-  .n_detected / .ncell
        .ratio
      }
      
    )
  ) |> 
  dplyr::select(-data) |>  
  dplyr::filter(.ratio > 0.5) ->
  flu2_variant_celltype


flu5_variant |> 
  dplyr::filter(sctype %in% inter_celltype) |> 
  dplyr::select(variant, af, sctype) |> 
  dplyr::group_by(variant) |> 
  tidyr::nest() |> 
  dplyr::ungroup() |> 
  dplyr::mutate(
    .ratio = purrr::map_dbl(
      .x = data,
      .f = function(.x) {
        .ncell <- nrow(.x)
        .x |> 
          dplyr::filter(!is.na(af)) |> 
          dplyr::filter(af != 0) |> 
          nrow() ->
          .n_detected
        
        .ratio <-  .n_detected / .ncell
        .ratio
      }
      
    )
  ) |> 
  dplyr::select(-data) |>  
  dplyr::filter(.ratio > 0.5) ->
  flu5_variant_celltype


ggvenn::ggvenn(
  data = list(
    Flu2 = unique(flu2_variant_celltype$variant),
    Flu5 = unique(flu5_variant_celltype$variant)
  ),
  show_percentage = FALSE,
  fill_color = ggsci::pal_aaas()(3),
  fill_alpha = 0.9,
  stroke_size = 0.5,
  set_name_size = 10,
  text_size = 8
)

ggplot() +  
  annotate(
    "text", x = 1, y = 1,
    size = 6,
    label = stringr::str_wrap(
      string = intersect(
        unique(flu2_variant_celltype$variant),
        unique(flu5_variant_celltype$variant)
      ) |> paste0(collapse = ", "),
      width = 30
    )
  ) + 
  theme_void()

intersect(
  unique(flu2_variant_celltype$variant),
  unique(flu5_variant_celltype$variant)
) |> paste0(collapse = ", ")

ggplot() +  
  annotate(
    "text", x = 1, y = 1,
    size = 6,
    color = ggsci::pal_aaas()(1),
    label = stringr::str_wrap(
      string = setdiff(
        unique(flu2_variant_celltype$variant),
        unique(flu5_variant_celltype$variant)
      ) |> paste0(collapse = ", "),
      width = 30
    )
  ) + 
  theme_void()

setdiff(
  unique(flu2_variant_celltype$variant),
  unique(flu5_variant_celltype$variant)
) |> paste0(collapse = ", ")

ggplot() +  
  annotate(
    "text", x = 1, y = 1,
    size = 6,
    color = ggsci::pal_aaas()(2)[2],
    label = stringr::str_wrap(
      string = setdiff(
        unique(flu5_variant_celltype$variant),
        unique(flu2_variant_celltype$variant)
      ) |> paste0(collapse = ", "),
      width = 30
    )
  ) + 
  theme_void()

setdiff(
  unique(flu5_variant_celltype$variant),
  unique(flu2_variant_celltype$variant)
) |> paste0(collapse = ", ")


# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------