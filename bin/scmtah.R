# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed Mar 22 15:04:30 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(ComplexHeatmap)
library(httr)

# src ---------------------------------------------------------------------
pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |> 
  dplyr::arrange(cancer_types)


# args --------------------------------------------------------------------

args <- commandArgs(TRUE)

barcode_cluster_file <- args[1]
cell_hetero_file <- args[2]
cell_coverage_file <- args[3]

cluster_hetero_file <- args[4]
cluster_coverage_file <- args[5]

cell_hetero_raw_file <- args[6]
# 
# barcode_cluster_file <- "/scr1/users/liuc9/tmp/mito/flu2-a/cromwell-executions/SCMTAH/b4545ece-8969-47e7-aea8-2dfeb2d1f872/call-cell_cluster_annotation/execution/barcode_cluster.tsv"
# cell_hetero_file <- "/scr1/users/liuc9/tmp/mito/flu2-a/cromwell-executions/SCMTAH/b4545ece-8969-47e7-aea8-2dfeb2d1f872/call-call_mt_variants/execution/cell/final/cell.cell_heteroplasmic_df.tsv.gz"
# cell_coverage_file <- "/scr1/users/liuc9/tmp/mito/flu2-a/cromwell-executions/SCMTAH/b4545ece-8969-47e7-aea8-2dfeb2d1f872/call-call_mt_variants/execution/cell/final/cell.coverage.txt.gz"
# 
# cluster_hetero_file <- "/scr1/users/liuc9/tmp/mito/flu2-a/cromwell-executions/SCMTAH/b4545ece-8969-47e7-aea8-2dfeb2d1f872/call-call_mt_variants/execution/cluster/final/cluster.cell_heteroplasmic_df.tsv.gz"
# cluster_coverage_file <- "/scr1/users/liuc9/tmp/mito/flu2-a/cromwell-executions/SCMTAH/b4545ece-8969-47e7-aea8-2dfeb2d1f872/call-call_mt_variants/execution/cluster/final/cluster.coverage.txt.gz"
# 
# cell_hetero_raw_file <- "/scr1/users/liuc9/tmp/mito/flu2-a/cromwell-executions/SCMTAH/b4545ece-8969-47e7-aea8-2dfeb2d1f872/call-call_mt_variants/execution/cell/final/cell.cell_heteroplasmic_df_raw.tsv.gz"


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
    dplyr::select(barcode, pos, depth) |>
    dplyr::arrange(pos) |>
    tidyr::pivot_wider(
      names_from = pos,
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

# cell cluster ------------------------------------------------------------

cluster_umap <- fn_load_cluster(
  .filename = barcode_cluster_file
)

# Cell allele -------------------------------------------------------------

cell_hetero <- fn_load_hetero(
  .filename = cell_hetero_file
)

cell_coverage <- fn_load_coverage(
  .filename = cell_coverage_file
)

cell_cluster_af <- fn_af(
  .cluster = cluster_umap, 
  .hetero = cell_hetero
)

cell_cluster_forplot <- fn_forplot(
  .af = cell_cluster_af, 
  .coverage = cell_coverage
)

ch_af_depth <- fn_heatmap(
  .forplot = cell_cluster_forplot
  )



{
  pdf(
    file = "cell_af_heatmap.pdf",
    width = 14, 
    height = 7
  )
  ComplexHeatmap::draw(object = ch_af_depth$ch_af)
  dev.off()

  pdf(
    file = "cell_depth_heatmap.pdf",
    width = 14, 
    height = 7
  )
  ComplexHeatmap::draw(object = ch_af_depth$ch_depth)
  dev.off()
}


# cluster allele-----------------------------------------------------------------


cluster_hetero <- fn_load_hetero(
  .filename = cluster_hetero_file
) |> 
  dplyr::mutate(cluster = barcode) |> 
  dplyr::mutate(cluster = factor(cluster)) |> 
  dplyr::left_join(
    cluster_umap |> 
      dplyr::mutate(cluster = celltype) |> 
      dplyr::mutate(cluster = factor(cluster)) |> 
      dplyr::select(cluster, celltype) |> 
      dplyr::distinct(),
    by = "cluster"
  ) |> 
  dplyr::select(-cluster) |> 
  dplyr::rename(cluster = celltype)

cluster_coverage <- fn_load_coverage(
  .filename = cluster_coverage_file
)


cluster_cluster_af <- 
  cluster_hetero |> tidyr::pivot_wider(
      names_from  = variant, 
      values_from = af
  )

cluster_cluster_forplot <- fn_forplot(
  .af = cluster_cluster_af, 
  .coverage = cluster_coverage
  )


cluster_ch_af_depth <- fn_heatmap(
  .forplot = cluster_cluster_forplot
  )

{
  pdf(
    file = "cluster_af_heatmap.pdf",
    width = 7, 
    height = 7
  )
  ComplexHeatmap::draw(object = cluster_ch_af_depth$ch_af)
  dev.off()
  
  pdf(
    file = "cluster_depth_heatmap.pdf",
    width = 7, 
    height = 7
  )
  ComplexHeatmap::draw(object = cluster_ch_af_depth$ch_depth)
  dev.off()
}

# Variant annotation ------------------------------------------------------


cluster_cluster_forplot$forplot |> 
  dplyr::filter(!is.na(depth)) |> 
  # dplyr::select(barcode, pos, variant) |>  
  dplyr::select(pos, variant) |> 
  dplyr::distinct() |> 
  dplyr::mutate(variant = gsub(
    pattern = "[0-9]*",
    replacement = "",
    x = variant
  )) |> 
  tidyr::separate(
    col = variant,
    into = c("ref", "var")
  ) |> 
  # dplyr::rename(sample = barcode) |> 
  dplyr::mutate(sample = "sample1") |> 
  dplyr::select(
    sample = sample, 
    pos = pos, 
    ref = ref,
    var = var
  ) ->
  cell_variants

readr::write_tsv(
  x = cell_variants,
  file = "cell_snvlist.tsv"
)

tryCatch(
  {
    cell_variant_response <- POST(
      "https://mitomap.org/mitomaster/websrvc.cgi",
      body = list(
        file = upload_file("cell_snvlist.tsv"),
        fileType = "snvlist",
        output = "detail"
      ),
      encode = "multipart"
    )
  },
  error = function(err) {
    print(paste("HTTP error:", err$message))
  },
  warning = function(w) {
    print(paste("Warning:", w$message))
  },
  finally = {
    print("Done.")
  }
)

cell_anno <- content(
  x = cell_variant_response,
  as = "text",
  encoding = "UTF-8"
) |> 
  data.table::fread(
    sep = "\t"
  )

readr::write_tsv(
  x = cell_anno,
  file = "cell_variant_annotation.tsv"
)

writexl::write_xlsx(
  x = cell_anno,
  path = "cell_variant_annotation.xlsx"
)


cell_anno |> 
  dplyr::mutate(
    variant = glue::glue("{tpos}{tnt}>{qnt}")
  ) |> 
  dplyr::select(
    variant, ntchange, calc_locus, patientphenotype,
    conservation, verbose_haplogroup
  ) |> 
  dplyr::mutate(
    calc_locus = gsub(
      pattern = "<br>.*",
      replace = "",
      x = calc_locus
    )
  ) |> 
  dplyr::mutate(
    conservation = gsub(
      pattern = "%",
      replacement = "",
      x = conservation
    )
  ) |> 
  dplyr::mutate(
    patientphenotype = stringr::str_wrap(
      stringr::str_to_sentence(string = patientphenotype),
      width = 10
    )
  ) |> 
  dplyr::mutate(conservation = as.numeric(conservation)) |> 
  dplyr::select(
    Ntchange = ntchange,
    Locus = calc_locus,
    Conservation = conservation,
    Haplogroup = verbose_haplogroup,
    Phenotype = patientphenotype
  ) ->
  variant_annotation
  


# Cluster cell allele -----------------------------------------------------

cell_hetero_raw <- fn_load_hetero(
  .filename = cell_hetero_raw_file
)

cell_raw_cluster_af <- cluster_umap |> 
  dplyr::left_join(cell_hetero_raw, by = "barcode") |>
  dplyr::rename(cluster = celltype) |> 
  dplyr::filter(variant %in% cluster_hetero$variant) |> 
  tidyr::pivot_wider(
    names_from = variant,
    values_from = af
  )
  
cell_raw_cluster_forplot <- fn_forplot(
  .af = cell_raw_cluster_af, 
  .coverage = cell_coverage
)

cell_raw_ch_af_depth <- fn_heatmap(
  .forplot = cell_raw_cluster_forplot, 
  .cell_variants = cell_cluster_forplot$forplot$variant,
  .variant_annotation = variant_annotation
  )

{
  pdf(
    file = "cluster_cell_af_heatmap.pdf",
    width = 14, 
    height = 7
  )
  ComplexHeatmap::draw(object = cell_raw_ch_af_depth$ch_af)
  dev.off()
  
  pdf(
    file = "cluster_cell_depth_heatmap.pdf",
    width = 14, 
    height = 7
  )
  ComplexHeatmap::draw(object = cell_raw_ch_af_depth$ch_depth)
  dev.off()
}




# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------

save.image(file = "scmtah.rda")
# load(file = "/scr1/users/liuc9/tmp/mito/flu2-a/cromwell-executions/SCMTAH/0138fcd0-c384-42c2-8704-6647767610d2/call-plot_scmtah/execution/scmtah.rda")
