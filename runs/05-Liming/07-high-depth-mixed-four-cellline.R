#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Mon Oct 23 15:58:26 2023
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
datadir <- "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth"

torundir <- file.path(
  datadir,
  "torun"
)

dir.create(torundir)

# body --------------------------------------------------------------------

tibble::tibble(
  path = list.files(
    file.path(datadir, "fastq"), 
    full.names = T
  )
) |> 
  dplyr::mutate(
    name = glue::glue(
      "{basename(path)}"
    )
  ) ->
  metadata

metadata |> 
  dplyr::mutate(
    a = purrr::map2(
      .x = path,
      .y = name,
      .f = \(.x, .y) {
        # .x <- metadata$path[[1]]
        # .y <- metadata$name[[1]]
        
        .ydir <- file.path(torundir, .y)
        dir.create(.ydir, showWarnings = F, recursive = T)
        
        .srrid <- .y
        
        conf <- list(
          "scMOCHA.output_id" = "{.srrid}" |> glue::glue(),
          "scMOCHA.fastqs" = "{.ydir}" |> glue::glue(),
          "scMOCHA.sample_id" = "{.srrid}" |> glue::glue(),
          "scMOCHA.transcriptome" = "/home/liuc9/data/refdata/mgatk_index/Human",
          "scMOCHA.rCRS" = "/home/liuc9/github/scMOCHA/fasta/rCRS.MT.fasta",
          "scMOCHA.output_dir" = "{.srrid}" |> glue::glue(),
          "scMOCHA.cellrefname" = "",
          "scMOCHA.celllevel" = "",
          "scMOCHA.memory" = "50 GB",
          "scMOCHA.boot_disk_size_gb" = "12",
          "scMOCHA.disk_space" = "50",
          "scMOCHA.cpu" = "10",
          "scMOCHA.scmocha_version" = "latest",
          "scMOCHA.docker" = "chunjiesamliu/scmocha",
          "scMOCHA.partition" = "defq",
          "scMOCHA.account" = "liuc9",
          "scMOCHA.IMAGE" = "/scr1/users/liuc9/sif/scmocha_latest.sif",
          "scMOCHA.perlscript" = "/home/liuc9/github/scMOCHA/bin/get_variants_info.pl",
          "scMOCHA.jar_path" = "/scr1/users/liuc9/tools/haplogrep3",
          "scMOCHA.sqlite_path" = "/mnt/isilon/xing_lab/liuc9/refdata/mitomaster/mitomap_sqlite_20230525.sqlite3"
        )
        
        conf_file <- file.path(
          .ydir,
          "{.srrid}.json" |> glue::glue()
        )
        
        jsonlite::write_json(
          x = conf,
          path = conf_file,
          auto_unbox = TRUE
        )
        
        .jsonfile <- file.path(
          .ydir, "{.srrid}.json" |> glue::glue()
        )
        .errfile <- file.path(
          .ydir, "{.srrid}.err" |> glue::glue()
        )
        .logfile <- file.path(
          .ydir, "{.srrid}.log" |> glue::glue()
        )
        runwdl_sh_file <- file.path(
          .ydir, "runwdl_{.srrid}.sh" |> glue::glue()
        )
        
        runwdl_cmd <- c(
          "#!/usr/bin/env bash",
          "# @AUTHOR: Chun-Jie Liu",
          "# @CONTACT: chunjie.sam.liu.at.gmail.com",
          "# @DATE: {lubridate::now()}" |> glue::glue(),
          "",
          "module load Java/15.0.1",
          "nohup java -Dconfig.file=/home/liuc9/github/scMOCHA/config/slurm.conf \\",
          "-jar /home/liuc9/tools/cromwell-78.jar \\",
          "run /home/liuc9/github/scMOCHA/scMOCHA-image.wdl \\",
          "-i {.jsonfile} 1>{.logfile} 2>{.errfile} &" |> glue::glue()
        )
        
        readr::write_lines(
          x = runwdl_cmd,
          file = runwdl_sh_file
        )
        
        tibble::tibble(
          path = list.files(.x, full.names = T) 
        ) |> 
          dplyr::mutate(
            targetname = basename(path)
          ) |> 
          dplyr::mutate(
            targetpath = file.path(.ydir, targetname)
          ) ->
          .xd
        
        .xd |> 
          dplyr::mutate(
            a = purrr::map2(
              .x = path,
              .y = targetpath,
              .f = \(.p, .tp) {
                tryCatch(
                  expr = {
                    file.symlink(
                      from = .p,
                      to = .tp
                    )
                  },
                  warning = \(w) {
                    print(w)
                  },
                  error = \(e) {
                    print(e)
                  }
                )
                
              }
            )
          )
        
        tibble::tibble(
          srrid = .y,
          srrdir = .ydir,
          
          scmocha_sh = runwdl_sh_file
        )
      }
    )
  ) |> 
  tidyr::unnest(cols = a) ->
  conf_scmocha

readr::write_csv(
  x = conf_scmocha,
  file = file.path(
    "/home/liuc9/github/scMOCHA/runs/05-Liming",
    "runfiles-high-depth.csv"
  )
)
# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------