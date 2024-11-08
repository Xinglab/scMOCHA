#!/usr/bin/env Rscript
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
# library(ComplexHeatmap)
suppressPackageStartupMessages(library(ComplexHeatmap))
library(httr)
library(GetoptLong)
library(logger)
ht_opt$message <- FALSE

# src ---------------------------------------------------------------------
pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
  dplyr::arrange(cancer_types)


# args --------------------------------------------------------------------


# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.
# 
# cell_meta_data_file <- "/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE157344/cromwell-executions/scMOCHABatch/87f7e7e9-4e27-491a-9125-19a78cddaf64/call-scMOCHA/shard-15/sub.scMOCHA/4d228407-529c-4b74-bc9e-f189ce7b2274/call-cell_cluster_annotation/execution/cell_meta_data.tsv"
# barcode_cluster_file <-"/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE157344/cromwell-executions/scMOCHABatch/87f7e7e9-4e27-491a-9125-19a78cddaf64/call-scMOCHA/shard-15/sub.scMOCHA/4d228407-529c-4b74-bc9e-f189ce7b2274/call-plot_scMOCHA/inputs/722149735/barcode_cluster.tsv"
# cell_hetero_file <-"/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE157344/cromwell-executions/scMOCHABatch/87f7e7e9-4e27-491a-9125-19a78cddaf64/call-scMOCHA/shard-15/sub.scMOCHA/4d228407-529c-4b74-bc9e-f189ce7b2274/call-plot_scMOCHA/inputs/2020544951/cell.cell_heteroplasmic_df.tsv.gz"
# cell_coverage_file <-"/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE157344/cromwell-executions/scMOCHABatch/87f7e7e9-4e27-491a-9125-19a78cddaf64/call-scMOCHA/shard-15/sub.scMOCHA/4d228407-529c-4b74-bc9e-f189ce7b2274/call-plot_scMOCHA/inputs/2020544951/cell.coverage.txt.gz"
# cluster_hetero_file <-"/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE157344/cromwell-executions/scMOCHABatch/87f7e7e9-4e27-491a-9125-19a78cddaf64/call-scMOCHA/shard-15/sub.scMOCHA/4d228407-529c-4b74-bc9e-f189ce7b2274/call-plot_scMOCHA/inputs/613722931/cluster.cell_heteroplasmic_df.tsv.gz"
# cluster_coverage_file <-"/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE157344/cromwell-executions/scMOCHABatch/87f7e7e9-4e27-491a-9125-19a78cddaf64/call-scMOCHA/shard-15/sub.scMOCHA/4d228407-529c-4b74-bc9e-f189ce7b2274/call-plot_scMOCHA/inputs/613722931/cluster.coverage.txt.gz"
# cell_hetero_raw_file <-"/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE157344/cromwell-executions/scMOCHABatch/87f7e7e9-4e27-491a-9125-19a78cddaf64/call-scMOCHA/shard-15/sub.scMOCHA/4d228407-529c-4b74-bc9e-f189ce7b2274/call-plot_scMOCHA/inputs/2020544951/cell.cell_heteroplasmic_df_raw.tsv.gz"
# perlscript <- "/home/liuc9/github/scMOCHA/bin/get_variants_info.pl"
# jar_path <- "/scr1/users/liuc9/tools/haplogrep3"
# sqlite_path <- "/mnt/isilon/xing_lab/liuc9/refdata/mitomaster/mitomap_sqlite_20230525.sqlite3"


conda_root <- "/home/liuc9/tools/anaconda3"
conda_env <- "scmocha"
verbose <- FALSE

spec <- "
Usage: Rscript scMOCHA.R [options]

Options:
<cell_meta_data_file|meta=s> cell_meta_data.tsv
<barcode_cluster_file=s> barcode_cluster.tsv
<cell_hetero_file|ceh=s> cell.cell_heteroplasmic_df.tsv.gz
<cell_coverage_file|cec=s> cell.coverage.txt.gz
<cluster_hetero_file|clh=s> cluster.cell_heteroplasmic_df.tsv.gz
<cluster_coverage_file|clc=s> cluster.coverage.txt.gz
<cell_hetero_raw_file|chr=s> cell.cell_heteroplasmic_df_raw.tsv.gz
<perlscript=s> /home/liuc9/github/scMOCHA/bin/get_variants_info.pl
<jar_path=s> /scr1/users/liuc9/tools/haplogrep3
<sqlite_path=s> /mnt/isilon/xing_lab/liuc9/refdata/mitomaster/mitomap_sqlite_20230525.sqlite3
<conda_root=s> /home/liuc9/tools/anaconda3
<conda_env=s> scmocha
<verbose!> Print messages
"

GetoptLong.options(help_style = "two-column")
GetoptLong(spec, template_control = list(opt_width = 50))



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
    ) |>
    dplyr::filter(af > 0.05) # filter variants which AF < 0.05
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

fn_load_meta <- function(.filename) {
  data.table::fread(
    input = .filename,
    sep = "\t"
  ) |> 
    dplyr::rename(
      barcode = cellbarcode
    ) |> 
    dplyr::select(-orig.ident)
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

fn_forplot <- function(.af, .coverage, .meta) {
  # print(.meta)
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
  
  .coverage |> 
    dplyr::group_by(barcode) |> 
    dplyr::summarise(sum_depth = sum(depth, na.rm = TRUE)) ->
    .coverage_cell

  list(
    rank = .rank,
    forplot = .forplot,
    meta = .meta,
    coverage_cell = .coverage_cell
  )
}


fn_heatmap <- function(.forplot, .cell_variants = NULL, .variant_annotation = NULL) {
  pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
    dplyr::arrange(cancer_types)
  
  .forplot$forplot |>
    dplyr::select(barcode, variant, af) |>
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

  .gcol <- if (is.null(.cell_variants)) {
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
    dplyr::left_join(
      .forplot$meta |> 
        dplyr::select(barcode, `MT%` = percent.mt),
      by = "barcode"
    ) |> 
    dplyr::left_join(
      .forplot$coverage_cell |> 
        dplyr::mutate(sum_depth = log10(sum_depth + 1)) |> 
        dplyr::rename(`log10(Total reads)` = sum_depth),
      by = "barcode"
    ) |> 
    tibble::column_to_rownames(var = "barcode") |>
    dplyr::rename(Cluster = cluster) ->
  .af_cluster


  col_clusters <- levels(.af_cluster$Cluster)
  col_colors <- pcc$color[1:length(levels(.af_cluster$Cluster))]

  names(col_colors) <- col_clusters

  chm_top <- ComplexHeatmap::HeatmapAnnotation(
    df = .af_cluster,
    # gap = unit(c(2, 2), "mm"),
    col = list(
      Cluster = col_colors,
      `MT%` = circlize::colorRamp2(
        breaks = c(2, 7, 10),
        colors = c("gold", "red", "black"),
        # colors =  c("#440154FF", "#FDE725FF"),
        space = "RGB"
      ),
      `log10(Total reads)` = circlize::colorRamp2(
        breaks =quantile(.af_cluster$`log10(Total reads)`, c(0.15, 0.75, 0.9), na.rm = T),
        colors = c("gold", "red", "blue"),
        space = "RGB"
      )
    ),
    which = "column"
  )

  ch_af <- if (!is.null(.variant_annotation)) {
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

metadata <- fn_load_meta(
  .filename = cell_meta_data_file
)
log_info("load metadata")
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
  .coverage = cell_coverage,
  .meta = metadata
)

log_info("fn_heatmap")
# print(cell_cluster_forplot)

ch_af_depth <- fn_heatmap(
  .forplot = cell_cluster_forplot
)



{
  pdf(
    file = "cell_af_heatmap.pdf",
    width = 14,
    height = 15
  )
  ComplexHeatmap::draw(object = ch_af_depth$ch_af)
  dev.off()

  pdf(
    file = "cell_depth_heatmap.pdf",
    width = 14,
    height = 15
  )
  ComplexHeatmap::draw(object = ch_af_depth$ch_depth)
  dev.off()
  log_success("save image")
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
  .coverage = cluster_coverage,
  .meta = metadata
)


cluster_ch_af_depth <- fn_heatmap(
  .forplot = cluster_cluster_forplot
)

{
  pdf(
    file = "cluster_af_heatmap.pdf",
    width = 7,
    height = 15
  )
  ComplexHeatmap::draw(object = cluster_ch_af_depth$ch_af)
  dev.off()

  pdf(
    file = "cluster_depth_heatmap.pdf",
    width = 7,
    height = 15
  )
  ComplexHeatmap::draw(object = cluster_ch_af_depth$ch_depth)
  dev.off()
}



# Cluster cell allele -----------------------------------------------------

cell_hetero_raw <- fn_load_hetero(
  .filename = cell_hetero_raw_file
) |>
  dplyr::filter(variant %in% cluster_hetero$variant)

cell_raw_cluster_af <- cluster_umap |>
  dplyr::left_join(cell_hetero_raw, by = "barcode") |>
  dplyr::rename(cluster = celltype) |>
  tidyr::pivot_wider(
    names_from = variant,
    values_from = af
  )

cell_raw_cluster_forplot <- fn_forplot(
  .af = cell_raw_cluster_af,
  .coverage = cell_coverage,
  .meta = metadata
)


# Variant annotation ------------------------------------------------------


cell_raw_cluster_forplot$forplot |>
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
  ) |>
  dplyr::mutate(
    v = glue::glue("{pos}{ref}>{var}")
  ) |>
  dplyr::select(sample, v) |>
  tibble::rowid_to_column() |>
  tidyr::pivot_wider(
    names_from = rowid,
    values_from = v
  ) ->
cell_variants

readr::write_delim(
  x = cell_variants,
  file = "cell_snvlist.tsv",
  delim = " ",
  col_names = F
)
# readr::write_tsv(
#   x = cell_variants,
#   file = "cell_snvlist.tsv"
# )

# fn_http_request <- function() {
#   cell_variant_response <- tryCatch(
#     {
#       POST(
#         "https://mitomap.org/mitomaster/websrvc.cgi",
#         body = list(
#           file = upload_file("cell_snvlist.tsv"),
#           fileType = "snvlist",
#           output = "detail"
#         ),
#         encode = "multipart"
#       )
#     },
#     error = function(err) {
#       print(paste("HTTP error:", err$message))
#       # "error"
#     },
#     warning = function(w) {
#       print(paste("Warning:", w$message))
#     },
#     finally = {
#       print("Done.")
#     }
#   )

#   status <- tryCatch(
#     expr = {
#       httr::status_code(cell_variant_response)
#     },
#     error = function(err) {
#       0
#     }
#   )



#   variant_annotation <- if(status == 200) {
#     cell_anno <- content(
#       x = cell_variant_response,
#       as = "text",
#       encoding = "UTF-8"
#     ) |>
#       data.table::fread(
#         sep = "\t"
#       )

#     readr::write_tsv(
#       x = cell_anno,
#       file = "cell_variant_annotation.tsv"
#     )

#     writexl::write_xlsx(
#       x = cell_anno,
#       path = "cell_variant_annotation.xlsx"
#     )


#     cell_anno |>
#       dplyr::mutate(
#         variant = glue::glue("{tpos}{tnt}>{qnt}")
#       ) |>
#       dplyr::select(
#         variant, ntchange, calc_locus, patientphenotype,
#         conservation, verbose_haplogroup
#       ) |>
#       dplyr::mutate(
#         calc_locus = gsub(
#           pattern = "<br>.*",
#           replace = "",
#           x = calc_locus
#         )
#       ) |>
#       dplyr::mutate(
#         conservation = gsub(
#           pattern = "%",
#           replacement = "",
#           x = conservation
#         )
#       ) |>
#       dplyr::mutate(
#         patientphenotype = stringr::str_wrap(
#           stringr::str_to_sentence(string = patientphenotype),
#           width = 30
#         )
#       ) |>
#       dplyr::mutate(conservation = as.numeric(conservation)) |>
#       dplyr::select(
#         Ntchange = ntchange,
#         Locus = calc_locus,
#         Conservation = conservation,
#         Haplogroup = verbose_haplogroup,
#         Phenotype = patientphenotype
#       )
#   } else {NULL}
# }

cmd <- "source {conda_root}/etc/profile.d/conda.sh; conda activate {conda_env}; perl {perlscript} {file.path(jar_path, 'haplogrep3.jar')} {sqlite_path} cell_snvlist.tsv > cell_variant_annotation.tsv" |> glue::glue()
# cmd <- "~/tools/anaconda3/envs/scmocha/bin/perl {perlscript} {file.path(jar_path, 'haplogrep3.jar')} {sqlite_path} cell_snvlist.tsv > cell_variant_annotation.tsv" |> glue::glue()
message(cmd)
system(command = cmd)

variant_annotation <- if (file.exists("cell_variant_annotation.tsv")) {
  cell_anno <- readr::read_tsv("cell_variant_annotation.tsv")
  writexl::write_xlsx(
    x = cell_anno,
    path = "cell_variant_annotation.xlsx"
  )


  cell_anno |>
    dplyr::mutate(
      variant = glue::glue("{Position}{Ref}>{Alt}")
    ) |>
    dplyr::mutate(
      Status = ifelse(
        !is.na(Status),
        "Reported",
        Status
      )
    ) |>
    dplyr::select(
      variant, ntchange,
      calc_locus = Locus,
      Haplogroup,
      Verbose_haplogroup,
      Disease,
      Status,
      Conservation,
      mito_freq = `Mitomap Frequency`,
      gnomad_freq = `Gnomad Frequency`
    ) |>
    dplyr::mutate(
      calc_locus = gsub(
        pattern = "<br>.*",
        replace = "",
        x = calc_locus
      )
    ) |>
    dplyr::mutate(
      Conservation = gsub(
        pattern = "%",
        replacement = "",
        x = Conservation
      )
    ) |>
    dplyr::mutate(
      Disease = stringr::str_wrap(
        stringr::str_to_sentence(string = Disease),
        width = 30
      )
    ) |>
    dplyr::mutate(Conservation = as.numeric(Conservation)) |>
    dplyr::mutate(
      mito_ref = mito_freq / 100,
      gnomad_freq = gnomad_freq / 100
    ) |>
    dplyr::select(
      Ntchange = ntchange,
      Locus = calc_locus,
      Haplogroup = Verbose_haplogroup,
      Disease = Disease,
      Status,
      Conservation,
      `Mitomap freq` = mito_freq,
      `Gnomad freq` = gnomad_freq
    )
} else {
  NULL
}


cell_raw_ch_af_depth <- fn_heatmap(
  .forplot = cell_raw_cluster_forplot,
  .cell_variants = cell_cluster_forplot$forplot$variant,
  .variant_annotation = variant_annotation
)

{
  pdf(
    file = "cluster_cell_af_heatmap.pdf",
    width = 20,
    height = 15
  )
  ComplexHeatmap::draw(object = cell_raw_ch_af_depth$ch_af)
  dev.off()

  pdf(
    file = "cluster_cell_depth_heatmap.pdf",
    width = 14,
    height = 15
  )
  ComplexHeatmap::draw(object = cell_raw_ch_af_depth$ch_depth)
  dev.off()
}




# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------

save.image(file = "scMOCHA.rda")
# load(file = "/scr1/users/liuc9/tmp/mito/flu2-a/cromwell-executions/SCMTAH/0138fcd0-c384-42c2-8704-6647767610d2/call-plot_scmtah/execution/scmtah.rda")
