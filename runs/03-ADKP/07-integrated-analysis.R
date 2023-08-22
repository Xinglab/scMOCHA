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

#library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------

outfiles <- readr::read_tsv(
  file = "/home/liuc9/github/scMOCHA/03-ADKP/output/outfiles.tsv"
)

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
        if(is.na(.x)) {return(NA)}
        
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
        if (.x == "FALSE") {return(NULL)}
        
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
        if(is.null(.x)) {return(NA_integer_)}
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
        if(is.null(.x)) {
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
    ratio = round(`number of cells after filtering`/ `estimated number of cells`,2)
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
    `# of cells`=`estimated number of cells`,
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
        if(is.na(.x)) {return(NULL)}
        
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
        
        if(is.null(.x)) {
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
    `Ratio of scMOCHA`= ratio1,
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
  dplyr::filter(!is.na(cellid)) |> 
  dplyr::count(azi_cluster, cluster) ->
  forplot


forplot |> 
  ggplot(aes(
    x = azi_cluster,
    y = cluster,
    fill = n,
    label = n
  )) +
  geom_label() +
  geom_tile() +
  scale_x_discrete(
    limit = paste("cluster", c(1:11, 13), sep = "_")
  ) +
  scale_y_discrete(
    limit = paste("cluster", c(1:11, 13), sep = "_")
  )

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "/scr1/users/liuc9/mitochondrial/realdata/03-ADKP/output/rda/07-integrated-analysis.rda")