# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed May  3 14:43:44 2023
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


# body --------------------------------------------------------------------

runfiles <- readr::read_tsv(
  file = "/home/liuc9/github/scMOCHA/runs/01-Sci_Immunol_32651212/runfiles.tsv"
)

runfiles |>
  dplyr::mutate(
    logfile = purrr::map_chr(
      .x = conf_json,
      .f = function(.x) {
        gsub(
          pattern = "json",
          replacement = "log",
          x = .x
        )
      }
    )
  ) |>
  dplyr::select(
    srrid, logfile
  ) |>
  dplyr::mutate(
    logfile_exists = file.exists(logfile)
  ) ->
  logfiles

logfiles |>
  dplyr::mutate(
    outfile = furrr::future_map_chr(
      .x = logfile,
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
            if(file.exists(.aa)) {
              .aa
            } else {
              FALSE
            }
          },
          error = function(err) {
            FALSE
          }
        )

      }
    )
  ) ->
  logfiles_outfile


logfiles_outfile |>
  dplyr::mutate(
    a = furrr::future_map2_chr(
      .x = logfile,
      .y = srrid,
      .f = function(.x, .y) {
        .xx <- readr::read_lines(file = .x)
        
        tryCatch(
          expr = {
            .xxx <- which(grepl(
              pattern = "SCMTAH.cluster_cell_af_heatmap",
              x = .xx
            ))[[2]]
            
            .a <- strsplit(.xx[.xxx], ":")[[1]][[2]]
            
            
            .aa <- gsub(
              pattern = " |\"|,",
              replacement = "",
              x = .a
            )
            .aaa <- file.path(dirname(.aa), "cell_snvlist.tsv")
            
            if(file.exists(.aaa)) {
              .aaa
            } else {
              FALSE
            }
          },
          error = function(err) {
            FALSE
          }
        )
        
      }
    )
  ) ->
  logfiles_outfile_a

logfiles_outfile_a |> 
  dplyr::mutate(
    b = purrr::map2(
      .x = srrid,
      .y = a,
      .f = function(.x, .y) {
        # .x <- logfiles_outfile_a$srrid[[1]]
        # .y <- logfiles_outfile_a$a[[1]]
        if(isFALSE(.y)) {
          return(NA)
        }
        
        .dir <- dirname(.y)
        .sfile <- file.path(
          .dir,
          "cell_snvlist.tsv"
        )
        .tfile <- file.path(
          "/home/liuc9/tmp/snvlist",
          glue::glue("{.x}_cell_snvlist.tsv")
        )
        
        file.copy(.sfile,.tfile)
        
        .sfile <- file.path(
          .dir,
          "cell_variant_annotation.xlsx"
        )
        .tfile <- file.path(
          "/home/liuc9/tmp/snvlist",
          glue::glue("{.x}_cell_variant_annotation.xlsx")
        )
        
        file.copy(.sfile,.tfile)
        
      }
    )
  )



logfiles_outfile |>
  dplyr::mutate(
    linkfile = purrr::map_chr(
      .x = outfile,
      .f = function(.x) {
        .xx <- basename(.x)

        if (.x == "FALSE") {
          return(NA)
        }

        .to <- "/home/liuc9/github/scMOCHA/01-Sci_Immunol_32651212/outputs/{.xx}" |> glue::glue()

        if(!file.exists(.to)) {
          file.symlink(
            from = .x,
            to = .to
          )
        }
        .to

      }
    )
  ) ->
  logfiles_outfile_linkfile

logfiles_outfile_linkfile |>
  dplyr::mutate(
    tardir = purrr::map_chr(
      .x = linkfile,
      .f = function(.x) {
        gsub(
          pattern = ".tar.gz",
          replacement = "",
          x = .x
        )
      }
    )
  ) |>
  dplyr::select(srrid, outfile, linkfile, tardir) ->
  outfiles

outfiles |>
  readr::write_tsv(
    file = "/scr1/users/liuc9/mitochondrial/realdata/01-Sci_Immunol_32651212/outputs/outfiles.tsv"
  )



# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------