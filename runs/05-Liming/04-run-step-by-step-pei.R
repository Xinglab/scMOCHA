#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Mon Oct 23 16:26:12 2023
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

runfiles <- data.table::fread(
  input =  file.path(
    "/home/liuc9/github/scMOCHA/runs/05-Liming",
    "runfiles_pei.csv"
  ),
  sep = ","
)


# body --------------------------------------------------------------------

runfiles$scmocha_sh |> 
  purrr::map_chr(
    .f = \(.x) {
      "bash {.x} &" |> glue::glue()
    }
  ) |> 
  readr::write_lines(
    file = file.path(
      "/home/liuc9/github/scMOCHA/05-Liming/cellline",
      "runwdl.sh"
    )
  )

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------