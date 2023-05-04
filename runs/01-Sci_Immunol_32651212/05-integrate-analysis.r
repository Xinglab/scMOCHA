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
    srrid, outfile, linkfile, tardir, Age, gender, 
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
  tidyr::unnest(cols = a) |> 
  dplyr::mutate(
    source_name = purrr::map2_chr(
      .x = source_name,
      .y = subject_status,
      .f = function(.x, .y) {
        if(is.na(.y)) {return(.x)}
        .y <- ifelse(
          .y == "Asymptomatic case of COVID-19 patient",
          yes = "mild COVID-19 patient",
          no = .y
        )
        .yy <- gsub(
          pattern = " COVID-19 patient",
          replacement = "",
          x = .y
        )
        glue::glue("{.x}({.yy})")
      }
    )
  ) ->
  metadata

metadata |> 
  readr::write_csv(
    file = "/home/liuc9/github/scRNAseq-MitoVariant/01-Sci_Immunol_32651212/outputs/metadata.csv"
  )


# body --------------------------------------------------------------------

metadata |> 
  dplyr::mutate(
    anno = purrr::map(
      .x = outfile,
      .f = function(.x) {
        if(.x == "FALSE") {return(NA)}
        
        .uuid <- dirname(dirname(dirname(.x)))
        
        .cva <- file.path(
          .uuid,
          "call-plot_scmtah/execution",
          "cell_variant_annotation.tsv"
        )
        
        readr::read_tsv(.cva, show_col_types = FALSE)
      }
    )
  ) |> 
  dplyr::mutate(
    nmut = purrr::map_int(
      .x = anno,
      .f = function(.x) {
        if(all(is.na(.x))) {return(NA_integer_)}
        nrow(.x)
      }
    )
  ) |> 
  dplyr::mutate(
    haplogroup = purrr::map(
      .x = anno,
      .f = function(.x) {
        if(all(is.na(.x))) {
          return(
            tibble::tibble(
              haplogroup = NA_character_,
              verbose_haplogroup = NA_character_
            )
          )
          
        }
        .x |> 
          dplyr::select(haplogroup, verbose_haplogroup) |> 
          dplyr::distinct()
      }
    )
  ) |> 
  tidyr::unnest(cols = haplogroup) ->
  metadata_anno


metadata_anno |> 
  dplyr::mutate(
    pass = ifelse(
      test = is.na(tardir),
      yes = "Fail",
      no = "Pass"
    )
  )|> 
  dplyr::select(srrid, Age, gender, source_name, subject_status, pass, `estimated number of cells`, `median UMI counts per cell`, `median genes per cell`, `number of cells after filtering`, haplogroup, verbose_haplogroup, nmut) |> 
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
    `Cell ratio` = ratio,
    `# of variants` = nmut,
    Haplogroup = haplogroup, 
    Haplogroup_v = verbose_haplogroup
  ) |> 
  dplyr::mutate(
    Status = gsub(
      pattern = " patient",
      replacement = "",
      x = Status
    )
  )  ->
  metadata_clean

metadata_clean |> 
  writexl::write_xlsx(
    path = "/scr1/users/liuc9/mitochondrial/realdata/01-Sci_Immunol_32651212/outputs/metadata_clean.xlsx"
  )


# metadata_anno -----------------------------------------------------------

cor.test(
  formula = ~ nmut + Age,
  data = metadata_anno |> 
    dplyr::filter(source_name == "Normal_PBMC")
)


cor.test(
  formula = ~ nmut + genderx,
  data = metadata_anno |> 
    dplyr::mutate(
      genderx = ifelse(
        gender == "male",
        1,
        0
      )
    )
)

metadata_anno |> 
  ggplot(
    aes(
      x = Age,
      y = nmut,
      color = source_name
    )
  ) +
  geom_point() +
  geom_smooth(method = "glm", se = FALSE) +
  # geom_line() +
  ggsci::scale_color_jama(
    name = "Source"
  ) +
  ggthemes::theme_base() +
  theme(
    axis.title = element_text(size = 32),
    legend.position = "top"
  ) +
  labs(
    x = "Age",
    y = "# of variants"
  ) ->
  age_cor_plot

ggsave(
  filename = "Age_nmut_cor.pdf",
  plo = age_cor_plot,
  device = "pdf",
  width = 9,
  height = 6,
  path = "/home/liuc9/github/scRNAseq-MitoVariant/01-Sci_Immunol_32651212/outputs"
)


metadata_anno |> 
  dplyr::filter(!is.na(linkfile)) |> 
  dplyr::mutate(
    variant = purrr::map_chr(
      .x = anno,
      .f = function(.x) {
        # .x |> 
        #   dplyr::mutate(variant = glue::glue("{tpos}{tnt}>{qnt}")) |> 
        #   dplyr::pull(variant)
        .x |> 
          dplyr::pull(verbose_haplogroup) |> 
          unique()
      }
    )
  ) |> 
  dplyr::select(srrid, source_name, variant) ->
  metadata_anno_v

metadata_anno_v |> 
  dplyr::filter(source_name == "Normal_PBMC") |> 
  dplyr::select(srrid, variant) |> 
  tibble::deframe() ->
  normal

metadata_anno_v |> 
  dplyr::filter(source_name == "Flu_PBMC") |> 
  dplyr::select(srrid, variant) |> 
  tibble::deframe() ->
  Flu_PBMC


# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------