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

pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |> 
  dplyr::arrange(cancer_types)


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
  )  |> 
  dplyr::slice(
    6:9, 1:5, 10:20
  ) ->
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


# cell ratio --------------------------------------------------------------

metadata_anno |> 
  dplyr::mutate(
    cellratio = purrr::map(
      .x = tardir,
      .f = function(.x) {
        if(is.na(.x)) {return(NULL)}
        .ratio <- readr::read_tsv(
          file = file.path(
            .x, "celltype_ratio.tsv"
          ),
          show_col_types = FALSE
        )
        .ratio
      }
    )
  ) ->
  metadata_anno_cellratio


metadata_anno_cellratio |> 
  dplyr::filter(!purrr::map_lgl(.x = cellratio, .f = is.null)) |> 
  dplyr::select(srrid, source_name, cellratio) |> 
  dplyr::mutate(color = dplyr::case_match(
    source_name,
    "Flu_PBMC" ~ggsci::pal_jama()(4)[[1]],
    "Normal_PBMC" ~ ggsci::pal_jama()(4)[[2]],
    "nCoV_PBMC(mild)" ~ ggsci::pal_jama()(4)[[3]],
    "nCoV_PBMC(severe)" ~ ggsci::pal_jama()(4)[[4]]
  )) |> 
  dplyr::slice(
    5:8, 1:4, 9:19
  ) |> 
  dplyr::arrange(dplyr::desc(dplyr::row_number())) ->
  for_ratio_plot

for_ratio_plot |> 
  tidyr::unnest(cellratio) |> 
  ggplot(aes(
    x = ratio,
    y = srrid,
    fill = celltype
  )) +
  geom_col() +
  ggsci::scale_fill_nejm(name = "Cell type") +
  scale_x_continuous(
    expand = expansion(mult = 0, add = 0),
    labels = scales::percent_format()
  ) +
  scale_y_discrete(
    limits = for_ratio_plot$srrid 
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size = 12, face = "bold"),
    axis.text.y = element_text(
      color = for_ratio_plot$color
    ),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Cell ratio") ->
  p_cellratio


ggsave(
  filename = "Cell_ratio.pdf",
  plo = p_cellratio,
  device = "pdf",
  width = 11,
  height = 5,
  path = "/home/liuc9/github/scRNAseq-MitoVariant/01-Sci_Immunol_32651212/outputs"
)


# Depth -------------------------------------------------------------------

metadata_anno_cellratio

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------