# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun Apr 30 18:55:54 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------


# body --------------------------------------------------------------------
runwdl <- readr::read_lines(
  file = "/home/liuc9/github/scMOCHA/02-Neuron_2022/runwdl.sh"
)

tibble::tibble(
  runwdl = runwdl
) |>
  dplyr::mutate(
    sh = gsub(
      pattern = "bash | &",
      replacement = "",
      x = runwdl
    )
  ) ->
  runwdl_sh


runwdl_sh |>
  dplyr::mutate(
    job = purrr::map_chr(
      .x = sh,
      .f = function(.x) {

        .srr <- gsub(
          pattern = "runwdl_|.sh",
          replacement = "",
          x = .x
        )
        glue::glue("{.srr}.log")
      }
    )
  ) ->
  runwdl_sh_log

runwdl_sh_log |>
  dplyr::mutate(
    tf = purrr::map_lgl(
      .x = job,
      .f = function(.x) {
        .xx <- readr::read_lines(file = .x)

        tryCatch(
          expr = {
            .xxx <- which(grepl(
              pattern = "SCMTAH.output_dir_tar_gz",
              x = .xx
            ))[[2]]

            .a <- strsplit(.xx[.xxx], ":")[[1]][[2]]


            .aa <- gsub(
              pattern = " |\"|,",
              replacement = "",
              x = .a
            )
            file.exists(.aa)
          },
          error = function(err) {
            FALSE
          }
        )

      }
    )
  ) ->
  runwdl_sh_log_tf

runwdl_sh_log_tf |>
  dplyr::filter(!tf) |>
  dplyr::select(runwdl) ->
  torun

readr::write_lines(
  x = torun$runwdl,
  file = "/home/liuc9/github/scMOCHA/02-Neuron_2022/new_runwdl.sh"
)

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------