#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Jan 25 15:18:47 2024
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
library(data.table)
#library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------


# body --------------------------------------------------------------------

m <- data.table::fread("/home/liuc9/github/scMOCHA/05-Liming/scmocha-celline/cromwell-executions/scMOCHA/c8470e56-2211-4779-8f58-ee9fc7ddc110/call-call_mt_variants/execution/cluster/final/cluster.coverage.txt.gz")
mm <- data.table::fread("/home/liuc9/github/scMOCHA/05-Liming/scmocha-celline/cromwell-executions/scMOCHA/c8470e56-2211-4779-8f58-ee9fc7ddc110/call-call_mt_variants/execution/cluster/final/cluster.cell_heteroplasmic_df.tsv.gz")


c(0, 1, 2, 3) |> 
  purrr::map(.f = \(.x) {
    data.table::fread(
      input = "/scr1/users/liuc9/mitochondrial/realdata/05-Liming/cluster{.x}.txt" |> glue::glue()
    )
  }) |> 
  dplyr::bind_rows() |> 
  dplyr::mutate(a = purrr::map(
    .x = Variant,
    .f = \(.x) {
      .s <- strsplit(.x, split = ">")[[1]]
      
      .n <- stringr::str_extract(.s[[1]], "\\d+") |> as.integer()
      .ref <- stringr::str_extract(.s[[1]], "[A-Za-z]+")
      
      .v <- .s[2]
      
      tibble::tibble(
        pos = .n,
        ref = .ref,
        variant = .v
      )
    }
  )) |> 
  tidyr::unnest(cols = a) ->
  shiping_variant

rcrs <- readr::read_lines(
  file = "/home/liuc9/github/scMOCHA/fasta/rCRS.MT.fasta",
  skip = 1
) |> 
  paste0(collapse = "") 
strsplit(rcrs, "")[[1]] |> 
  tibble::enframe(name = "pos", "ref") ->
  rcrs_ref

bases <- c("A", "G", "C", "T")
names(bases) <- bases

base_list <- as.list(bases)

base_list |> 
  purrr::map(
    .f = \(.x) {
      data.table::fread(
        input = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-celline/cromwell-executions/scMOCHA/c8470e56-2211-4779-8f58-ee9fc7ddc110/call-call_mt_variants/execution/cluster/final/cluster.{.x}.txt.gz" |> glue::glue(),
        col.names = c("pos", "cluster", "f", "b")
      ) |> 
        dplyr::mutate(totalcount = f + b) |> 
        dplyr::mutate(variant = .x) |> 
        dplyr::select(pos, cluster, totalcount, variant) |>  
        tidyr::spread(key = cluster, value = totalcount)  ->
        .d
      
      rcrs_ref |> 
        dplyr::left_join(.d, by = "pos") |> 
        tidyr::replace_na(
          replace = list(
            cluster_0 = 0,
            cluster_1 = 0,
            cluster_2 = 0,
            cluster_3 = 0
          )
        )
    }
  ) |> 
  dplyr::bind_rows() ->
  base_list_load



shiping_variant |> 
  dplyr::group_by(cluster) |> 
  tidyr::nest() |> 
  dplyr::ungroup() |> 
  dplyr::mutate(cluster = gsub("cluster", "cluster_", cluster)) |> 
  dplyr::mutate(
    p = purrr::map2(
      .x = data,
      .y = cluster,
      .f = \(.x, .y) {
        .x |> 
          dplyr::mutate(aref = glue::glue("{pos}{ref}")) |> 
          dplyr::arrange(pos) |> 
          dplyr::mutate(aref = factor(aref, levels = aref)) |> 
          dplyr::arrange(variant) |> 
          dplyr::mutate(variant = factor(variant)) -> 
          .xd
        .xd |> 
          ggplot(aes(x = aref, y = variant)) +
          geom_tile(aes(fill = Frequency)) +
          geom_label(aes(label = Frequency)) +
          theme_bw() +
          theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, face = "bold"),
            axis.text.y = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 16, face = "bold")
          ) +
          labs(
            x = "Reference",
            y = "Variant"
          ) ->
          p_shiping
        
        base_list_load |> 
          dplyr::select(pos, ref, variant, Count = .y) |> 
          dplyr::filter(pos %in% .xd$pos) |> 
          dplyr::filter(!is.na(variant)) |> 
          dplyr::mutate(aref = glue::glue("{pos}{ref}")) |> 
          dplyr::arrange(pos) |> 
          dplyr::mutate(aref = factor(aref, levels = unique(aref))) |> 
          dplyr::arrange(variant) |> 
          dplyr::mutate(variant = factor(variant)) ->
          .xxd
        
        .xxd |> 
          ggplot(aes(x = aref, y = variant)) +
          geom_tile(aes(fill = Count)) +
          geom_label(aes(label = Count)) +
          theme_bw() +
          theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, face = "bold"),
            axis.text.y = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            axis.title.x = element_blank()
          ) +
          labs(
            x = "Reference",
            y = "Variant"
          )  ->
          p_cj
        
        p_cj / p_shiping
        
      }
    )
  ) ->
  shiping_variant_cluster


shiping_variant_cluster |> 
  dplyr::mutate(
    a = purrr::map2(
      .x = cluster,
      .y = p,
      .f = \(.x, .y) {
        ggsave(
          plot = .y,
          filename = "freq_count_{.x}.pdf" |> glue::glue(),
          device = "pdf",
          path = "/home/liuc9/github/scMOCHA/05-Liming/benchmark",
          width = 12,
          height = 7
        )
      }
    )
  )
# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------