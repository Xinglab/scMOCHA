# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Apr 27 16:47:33 2023
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

datadir <- "/home/liuc9/scratch/mitochondrial/realdata/02-Neuron_2022"

runfiles <- readr::read_tsv(
  file = file.path(
    datadir, "runfiles.tsv"
  )
)

# step1 -------------------------------------------------------------------

runfiles$dump_slrm |>
  purrr::map_chr(
    .f = function(.x) {
      command <- "/cm/shared/apps/slurm/current/bin/sbatch {.x}" |> glue::glue()
      # system(command = command)
    }
  ) |> 
  readr::write_lines(
    file = file.path(
      datadir, "sra_dump.sh"
    )
  )

# step2 -------------------------------------------------------------------



runfiles$scmtah_sh |>
  purrr::map_chr(
    .f = function(.x) {
      "bash {.x} &" |> glue::glue()
      # system(command = cmd)
    }
  ) |> 
  readr::write_lines(
    file = file.path(
      datadir, "runwdl.sh"
    )
  )
  


# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------