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
    # dplyr::filter(source_name == "Normal_PBMC"),
    dplyr::filter(source_name == "nCoV_PBMC(severe)"),
  method = "spearman"
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
  # ggthemes::theme_base() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16, face = "bold", colour = "black"),
  ) +
  labs(
    x = "Age",
    y = "# of variants"
  ) ->
  age_cor_plot;age_cor_plot

ggsave(
  filename = "Age_nmut_cor.pdf",
  plo = age_cor_plot,
  device = "pdf",
  width = 9,
  height = 6,
  path = "/home/liuc9/github/scRNAseq-MitoVariant/01-Sci_Immunol_32651212/outputs"
)

metadata_anno |> 
  ggplot(
    aes(
      x = `number of cells after filtering`,
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
  # ggthemes::theme_base() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16, face = "bold", colour = "black")
  ) +
  labs(
    x = "# of cells",
    y = "# of variants"
  ) ->
  ncells_cor_plot

ggsave(
  filename = "Ncells_nmut_cor.pdf",
  plo = ncells_cor_plot,
  device = "pdf",
  width = 9,
  height = 6,
  path = "/home/liuc9/github/scRNAseq-MitoVariant/01-Sci_Immunol_32651212/outputs"
)




cor.test(
  formula = ~ nmut + Age,
  data = metadata_anno |> 
    # dplyr::filter(source_name == "Normal_PBMC"),
    dplyr::filter(source_name == "nCoV_PBMC(severe)"),
  method = "spearman"
)



wilcox.test(
  formula = nmut ~ gender ,
  data = metadata_anno |> 
    dplyr::filter(!is.na(linkfile)) |> 
    dplyr::mutate(gender = factor(gender))
)

t.test(
  formula = nmut ~ gender ,
  data = metadata_anno |> 
    dplyr::filter(!is.na(linkfile)) |> 
    dplyr::mutate(gender = factor(gender))
) |> broom::tidy()


metadata_anno |> 
  ggplot(aes(
    x = gender,
    y = nmut,
  )) +
  # geom_violin() +
  geom_boxplot(
    aes(fill = gender),
    width = 0.5,
    show.legend = FALSE
  ) +
  geom_point(position = position_jitter(width = 0.3)) +
  theme_bw() +
  ggsci::scale_fill_aaas() +
  scale_x_discrete(
    limits = c("male", "female"),
    labels = c("Male", "Female")
  ) +
  theme(
    axis.title = element_text(size = 16, face = "bold", colour = "black"),
    axis.text.x = element_text(size = 14, face = "bold", colour = "black"),
    axis.title.x = element_blank()
  ) +
  labs(
    x = "Gender",
    y = "# of variants"
  ) ->
  gender_cor_plot


(age_cor_plot / ncells_cor_plot | gender_cor_plot) +
  plot_layout(
    width = c(3, 1),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(legend.position = "bottom") ->
  age_ncells_p;age_ncells_p

ggsave(
  filename = "Ncells_age_nmut_cor.pdf",
  plo = age_ncells_p,
  device = "pdf",
  width = 10,
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

gtf_gene_df <- 
  readr::read_rds(
    file = "/home/liuc9/github/scRNAseq-MitoVariant/fasta/mt_exons.df.rds.gz"
  )

library(ggtranscript)
gtf_gene_df %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = gene_name
  )) +
  geom_range( aes(fill = transcript_biotype)) +
  geom_intron(
    data = to_intron(gtf_gene_df, "transcript_name"),
    aes(strand = strand)
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.03)),
    limits = c(1, 17000),
    breaks = seq(1000, 17000, 1000),
    labels = seq(1000, 17000, 1000)
  ) +
  # scale_fill_brewer(palette = "Set3")
  ggsci::scale_fill_jama(
    name = "Biotype",
    labels = c("MT rRNA", "MT tRNA", "Protein coding")
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    axis.line = element_line(color = "black"),
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.title.x = element_blank(),
    # axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    x = "Position"
  ) ->
  p_mt_chrom

metadata_anno_cellratio |> 
  dplyr::mutate(
    depth = purrr::map(
      .x = tardir,
      .f = function(.x) {
        if(is.na(.x)) {return(NULL)}
        data.table::fread(
          input = file.path(
            .x, "possorted_genome_bam.MT.depth"
          ),
          col.names =  c("chr", "pos", "depth")
        )
        
      }
    )
  ) ->
  metadata_anno_depth

metadata_anno_depth |> 
  dplyr::filter(!purrr::map_lgl(.x = depth, .f = is.null)) |> 
  dplyr::select(srrid, source_name, depth) |> 
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
  dplyr::mutate(srrid = factor(srrid, levels = srrid)) ->
  for_depth_plot

for_depth_plot |> 
  tidyr::unnest(cols = depth) |>
  ggplot(aes(x=pos, y = depth, fill = srrid)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.03)),
    limits = c(1, 17000),
  ) +
  scale_y_continuous(
    expand = c(0.01, 0),
    # label = scales::label_number(scale = 1e-5, suffix = "x10^5")
  ) +
  scale_fill_manual(
    name = "Sample",
    values = for_depth_plot$color,
    guide = guide_legend(nrow = 3)
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    # axis.title.y = element_blank(),
    axis.line.x.bottom = element_line(color = "black"),
    # strip.background = element_rect(fill = NA, colour = "black"),
    strip.background = element_blank(),
    # strip.text = element_text(
    #   color = "black",
    #   face = "bold",
    #   size = 8
    # ),
    # strip.text = element_text(
    #   color = for_depth_plot$color
    # ),
    strip.text = element_blank(),
    legend.position = "none"
  ) +
  facet_wrap(
    facets = ~srrid,
    ncol = 1,
    strip.position = "right"
  ) +
  labs(y = "Depth") ->
  p_mt_depth

ggsave(
  filename = "Sample_depth.pdf",
  plo = p_mt_depth,
  device = "pdf",
  width = 15,
  height = 10,
  path = "/home/liuc9/github/scRNAseq-MitoVariant/01-Sci_Immunol_32651212/outputs"
)

p_depth <- cowplot::plot_grid(
  plotlist = list(p_mt_depth, p_mt_chrom),
  ncol = 1,
  align = "v",
  rel_heights = c(0.7, 0.3)
)

ggsave(
  filename = "Sample_depth_merge.pdf",
  plo = p_depth,
  device = "pdf",
  width = 15,
  height = 15,
  path = "/home/liuc9/github/scRNAseq-MitoVariant/01-Sci_Immunol_32651212/outputs"
)



# Normal ------------------------------------------------------------------

metadata_anno_depth |> 
  readr::write_rds(
    file = "/home/liuc9/github/scRNAseq-MitoVariant/01-Sci_Immunol_32651212/outputs/metadata_anno_depth.rds"
  )


# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image("/home/liuc9/github/scRNAseq-MitoVariant/01-Sci_Immunol_32651212/outputs/05-integrate-analysis.rda")

# load(file = "/home/liuc9/github/scRNAseq-MitoVariant/01-Sci_Immunol_32651212/outputs/05-integrate-analysis.rda")
