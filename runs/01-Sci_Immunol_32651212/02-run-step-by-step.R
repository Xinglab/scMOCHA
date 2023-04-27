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

runfiles <- readr::read_tsv(
  file = "runs/01-Sci_Immunol_32651212/runfiles.tsv"
)

# step1 -------------------------------------------------------------------

runfiles$dump_slrm |> 
  purrr::map(
    .f = function(.x) {
      command <- "/cm/shared/apps/slurm/current/bin/sbatch {.x}" |> glue::glue()
      # print(cmd)
      system(command = command)
    }
  )

# step2 -------------------------------------------------------------------



# runfiles$scmtah_sh |> 
#   purrr::map(
#     .f = function(.x) {
#       cmd <- "bash {.x}" |> glue::glue()
#       system(command = cmd)
#     }
#   )


# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------