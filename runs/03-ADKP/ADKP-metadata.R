# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed May 24 13:58:48 2023
# @DESCRIPTION:


# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(httr)
library(synapser)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------


# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
synLogin("chunjie.sam.liu", "!Uu201012670")

# body --------------------------------------------------------------------
query <- synTableQuery("SELECT * FROM syn11346063.37")
query$filepath
csv_filepath <- "/scr1/users/liuc9/mitochondrial/realdata/03-ADKP/syn11346063.37.csv"
# file.copy(query$filepath, csv_filepath, overwrite = T)
file.exists(csv_filepath)


meta <- vroom::vroom(
  csv_filepath,
  delim = ","
) |> 
  data.table::as.data.table()

meta

meta |> 
  dplyr::select(
    id, name, study, dataType, assay, organ, tissue, species, sex, consortium, modelSystemName, treatmentType, specimenID, individualID, individualIdSource, specimenIdSource, resourceType, dataSubtype, metadataType, assayTarget, analysisType, cellType, nucleicAcidSource, fileFormat, group, projectId, libraryPrep
  )  ->
  meta_sel


meta_sel$grant |> table()

meta_sel |> 
  dplyr::filter(id == "syn21389259")

meta_sel |> 
  dplyr::filter(projectId == "syn2580853") |> 
  dplyr::filter(species == '["Human"]') |> 
  dplyr::filter(!is.na(individualID)) ->
  meta_sel_fil


meta_sel_fil |> 
  dplyr::select(
    study, dataType, assay, organ, tissue, sex, consortium,
    specimenID, individualID, libraryPrep
  ) |> 
  dplyr::distinct() ->
  meta_sel_fil_ind

meta_sel_fil_ind |> 
  dplyr::group_by(individualID) |> 
  dplyr::count() |> 
  dplyr::arrange(-n)

meta_sel_fil_ind |> 
  dplyr::mutate_all(
    .funs = function(.x) {
      gsub(
        pattern = '"|\\[|\\]',
        replacement = "",
        x = .x
      )
    }
  ) ->
  meta_sel_fil_ind_r

meta_sel_fil_ind_r


meta_sel_fil_ind_r |> 
  dplyr::count(
    consortium, study, dataType, organ
  ) |> 
  dplyr::mutate(
    organ_r = n / sum(n)
  ) |> 
  dplyr::group_by(dataType) |> 
  dplyr::mutate(
    dataType_r = sum(organ_r)
  ) |> 
  dplyr::ungroup() |> 
  dplyr::group_by(study) |> 
  dplyr::mutate(
    study_r = sum(organ_r)
  ) |> 
  dplyr::ungroup() |> 
  dplyr::group_by(consortium) |> 
  dplyr::mutate(
    consortium_r = sum(organ_r)
  ) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(
    consortium = glue::glue("{consortium} {round(consortium_r * 100, 2)}%"),
    study = glue::glue("{study} {round(study_r * 100, 2)}%"),
    dataType = glue::glue("{dataType} {round(dataType_r * 100, 2)}%"),
    organ = glue::glue("{organ} {round(organ_r * 100, 2)}%")
  ) |> 
  dplyr::select(1,2,3,4,5) |> 
  plotme::count_to_sunburst() ->
  psun

reticulate::py_run_string("import sys")

plotly::save_image(
  p = psun,
  file = file.path(
    "/home/liuc9/github/scMOCHA/03-ADKP/output",
    glue::glue("metadata-sunburst.pdf")
  ),
  width = 800,
  height = 800,
  device = "pdf"
)

htmlwidgets::saveWidget(
  widget = psun,
  file = file.path(
    "/home/liuc9/github/scMOCHA/03-ADKP/output",
    glue::glue("metadata-sunburst.html")
  )
)


# scRNA -------------------------------------------------------------------


meta_sel_fil_ind_r |> 
  dplyr::count(dataType, assay) |> 
  dplyr::arrange(-n)

meta_sel_fil_ind_r |> 
  dplyr::count(assay) |> 
  dplyr::arrange(assay)

meta_sel_fil_ind_r |> 
  dplyr::filter(grepl(
    pattern = "geneExpression",
    x = dataType
  )) |> 
  dplyr::count(assay) |> 
  dplyr::arrange(-n)

meta_sel_fil_ind_r |> 
  dplyr::filter(assay == "scrnaSeq") |> 
  dplyr::count(
    consortium, study
  ) 

meta_sel_fil_ind_r |> 
  dplyr::filter(assay == "scrnaSeq") |> 
  dplyr::count(
    consortium, study, organ
  ) |> 
  dplyr::mutate(
    organ_r = n / sum(n),
  ) |>  
  dplyr::group_by(study) |> 
  dplyr::mutate(
    study_r = sum(organ_r),
    study_r_n = sum(n)
  ) |> 
  dplyr::ungroup() |> 
  dplyr::group_by(consortium) |> 
  dplyr::mutate(
    consortium_r = sum(organ_r),
    consortium_r_n = sum(n)
  ) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(
    consortium = glue::glue("{consortium}\n{round(consortium_r * 100, 2)}%, n={consortium_r_n}"),
    study = glue::glue("{study}\n{round(study_r * 100, 2)}%, n={study_r_n}"),
    organ = glue::glue("{organ}\n{round(organ_r * 100, 2)}%, n={n}")
  ) |>  
  dplyr::select(1,2,3,4) |> 
  plotme::count_to_sunburst() ->
  psun_filter;psun_filter

reticulate::py_run_string("import sys")
plotly::save_image(
  p = psun_filter,
  file = file.path(
    "/home/liuc9/github/scMOCHA/03-ADKP/output",
    glue::glue("metadata-sunburst-filter.pdf")
  ),
  width = 800,
  height = 800,
  device = "pdf"
)

htmlwidgets::saveWidget(
  widget = psun_filter,
  file = file.path(
    "/home/liuc9/github/scMOCHA/03-ADKP/output",
    glue::glue("metadata-sunburst-filter.html")
  )
)



# Save selected data ------------------------------------------------------

meta_sel_fil |>
  # dplyr::mutate(
  #   assay = purrr::map_chr(
  #     .x = assay, 
  #     .f = function(.x) {
  #       gsub(
  #         pattern = '"|\\[|\\]',
  #         replacement = "",
  #         x = .x
  #       )
  #     }
  #   )
  # ) |> 
  dplyr::mutate_all(
    .funs = function(.x) {
      gsub(
        pattern = '"|\\[|\\]',
        replacement = "",
        x = .x
      )
    }
  ) |> 
  dplyr::filter(assay == "scrnaSeq") |> 
  dplyr::filter(
    fileFormat == "fastq"
  ) |> 
  dplyr::filter(
    nucleicAcidSource == "single cell"
  ) |> 
  readr::write_csv(
    file = "/home/liuc9/github/scMOCHA/03-ADKP/output/selected-ad-samples.csv"
  )


# footer ------------------------------------------------------------------


# save image --------------------------------------------------------------
