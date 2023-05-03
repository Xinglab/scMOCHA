# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed May  3 15:05:40 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------

outfiles <- readr::read_tsv(
  file = "/scr1/users/liuc9/mitochondrial/realdata/01-Sci_Immunol_32651212/outputs/outfiles.tsv"
)

srarun <- readr::read_csv(
  file = "/scr1/users/liuc9/mitochondrial/realdata/01-Sci_Immunol_32651212/SraRunTable.txt"
) |> 
  dplyr::rename(
    srrid = Run
  )

outfiles |> 
  dplyr::left_join(
    srarun, 
    by = "srrid"
  ) ->
  outfiles_sra


outfiles_sra |> 
  dplyr::select(
    srrid, tardir, Age, gender, 
    samplename = `Sample Name`,
    source_name,
    subject_group,
    subject_status
  ) |> 
  dplyr::arrange(source_name, subject_status) |> 
  dplyr::mutate(
    a = purrr::map(
      .x = tardir,
      .f = function(.x) {
        if(is.na(.x)) {return(NA)}
        
        readxl::read_excel(
          path = file.path(
            .x,
            "qc_cell_stats.xlsx"
          )
        )
      }
    )
  ) |> 
  tidyr::unnest(cols = a) ->
  metadata

metadata |> 
  readr::write_csv(
    file = "/home/liuc9/github/scRNAseq-MitoVariant/01-Sci_Immunol_32651212/outputs/metadata.csv"
  )

# body --------------------------------------------------------------------

metadata |> 
  dplyr::mutate(
    pass = ifelse(
      test = is.na(tardir),
      yes = "Fail",
      no = "Pass"
    )
  )|> 
  dplyr::select(srrid, Age, gender, source_name, subject_status, pass, `estimated number of cells`, `median UMI counts per cell`, `median genes per cell`, `number of cells after filtering`) |> 
  dplyr::arrange(
    source_name, 
    subject_status,
    pass,
    Age
  ) |> 
  dplyr::mutate(
    ratio = round(`number of cells after filtering`/ `estimated number of cells`,2)
  ) |> 
  dplyr::select(-pass) |> 
  dplyr::select(
    SRRID = srrid,
    Age,
    Gender = gender,
    Source = source_name,
    Status = subject_status,
    `Median UMI/cell` = `median UMI counts per cell`,
    `Median genes/cell` = `median genes per cell`,
    `# of cells`=`estimated number of cells`,
    `# cells after filter` = `number of cells after filtering`,
    `Cell ratio` = ratio
  ) |> 
  dplyr::mutate(
    Status = gsub(
      pattern = " patient",
      replacement = "",
      x = Status
    )
  ) ->
  metadata_clean

metadata_clean |> 
  writexl::write_xlsx(
    path = "/scr1/users/liuc9/mitochondrial/realdata/01-Sci_Immunol_32651212/outputs/metadata_clean.xlsx"
  )



metadata$tardir[[2]]

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------