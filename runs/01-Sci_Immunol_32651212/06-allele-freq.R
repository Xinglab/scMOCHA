# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu May  4 15:46:02 2023
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
metadata_anno_depth <- 
  readr::read_rds(
    file = "/home/liuc9/github/scRNAseq-MitoVariant/01-Sci_Immunol_32651212/outputs/metadata_anno_depth.rds"
  )


# body --------------------------------------------------------------------

metadata_anno_depth |> dplyr::glimpse()

metadata_anno_depth |> 
  dplyr::mutate(
    variant = purrr::map2(
      .x = anno,
      .y = tardir,
      .f = function(.x, .y) {
        if(is.na(.y)) {return(NULL)}
        .x |> 
          dplyr::mutate(
            variant = glue::glue("{tpos}{tnt}>{qnt}")
          ) |> 
          dplyr::select(variant)
      }
    )
  ) ->
  metadata_anno_depth_variant


metadata_anno_depth_variant |> 
  dplyr::mutate(color = dplyr::case_match(
    source_name,
    "Flu_PBMC" ~ggsci::pal_jama()(4)[[1]],
    "Normal_PBMC" ~ ggsci::pal_jama()(4)[[2]],
    "nCoV_PBMC(mild)" ~ ggsci::pal_jama()(4)[[3]],
    "nCoV_PBMC(severe)" ~ ggsci::pal_jama()(4)[[4]]
  )) |> 
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
    dplyr::select(-data) |> 
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
    path = "/home/liuc9/github/scRNAseq-MitoVariant/01-Sci_Immunol_32651212/outputs",
    width = 9, 
    height = 6,
    device = "pdf"
  )
  
  .p_up
}

metadata_anno_depth_variant$source_name |> 
  unique() |> 
  purrr::map(
    .f = fn_upset_plot
  ) ->
  p_ups

(p_ups[[2]] | p_ups[[1]]) / (p_ups[[3]] | p_ups[[4]]) +
  plot_annotation(tag_levels = "A") ->
  p_ups_together;p_ups_together
  
ggsave(
  plot = p_ups_together,
  filename = "upset-all.pdf" |> glue::glue(),
  path = "/home/liuc9/github/scRNAseq-MitoVariant/01-Sci_Immunol_32651212/outputs",
  width = 16, 
  height = 10,
  device = "pdf"
)


# for_variant |> 
#   dplyr::select(srrid, variant) |> 
#   tibble::deframe() |> 
#   purrr::reduce(.f = union) -> all_variants

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
  purrr::reduce(.f = intersect) -> common_variants

ggplot() +  
  annotate(
    "text", x = 1, y = 1,
    size = 6,
    color = "black", #ggsci::pal_aaas()(1),
    label = stringr::str_wrap(
      string = common_variants |> paste0(collapse = ", "),
      width = 30
    )
  ) + 
  theme_void()

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------

save.image(
  file = "/home/liuc9/github/scRNAseq-MitoVariant/01-Sci_Immunol_32651212/outputs/06-allele-freq.rda"
)
