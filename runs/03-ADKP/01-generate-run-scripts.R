# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Aug  3 17:07:54 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
#library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------

synid <- "syn12514624"
synid <- "syn21438358"
# load data ---------------------------------------------------------------

datadir <- "/scr1/users/liuc9/mitochondrial/realdata/03-ADKP"

metadata <- data.table::fread(
  input = file.path(
    datadir,
    "fastq",
    "SYNAPSE_METADATA_MANIFEST.tsv"
  ),
  sep = "\t"
)

metadata |> 
  dplyr::glimpse()

# body --------------------------------------------------------------------

metadata |> 
  dplyr::select(path, name, specimenID, individualID) |> 
  dplyr::arrange(specimenID, name)

torundir <- file.path(
  datadir,
  "torun"
)

metadata |> 
  dplyr::select(path, name, specimenID, diagnosis, individualID) |> 
  dplyr::group_by(diagnosis,specimenID, individualID) |> 
  tidyr::nest() |> 
  dplyr::ungroup() ->
  metadata_d

metadata_d |> 
  dplyr::mutate(
    a = purrr::map2(
      .x = data,
      .y = individualID,
      .f = \(.x, .y) {
        # .y <- metadata_d$individualID[[1]]
        # .x <- metadata_d$data[[1]]
        
        .ydir <- file.path(torundir, .y)
        print(.y)
        dir.create(.ydir, showWarnings = F, recursive = T)
        
        
        .srrid <- .y
        
        conf <- list(
          "scMOCHA.output_id" = "{.srrid}" |> glue::glue(),
          "scMOCHA.fastqs" = "{.ydir}" |> glue::glue(),
          "scMOCHA.sample_id" = "{.srrid}" |> glue::glue(),
          "scMOCHA.transcriptome" = "/home/liuc9/data/refdata/mgatk_index/Human",
          "scMOCHA.rCRS" = "/home/liuc9/github/scMOCHA/fasta/rCRS.MT.fasta",
          "scMOCHA.output_dir" = "{.srrid}" |> glue::glue(),
          "scMOCHA.cellrefname" = "mousecortexref",
          "scMOCHA.celllevel" = "subclass",
          "scMOCHA.memory" = "50 GB",
          "scMOCHA.boot_disk_size_gb" = "12",
          "scMOCHA.disk_space" = "50",
          "scMOCHA.cpu" = "10",
          "scMOCHA.scmocha_version" = "latest",
          "scMOCHA.docker" = "chunjiesamliu/scmocha",
          "scMOCHA.partition" = "defq",
          "scMOCHA.account" = "liuc9",
          "scMOCHA.IMAGE" = "/scr1/users/liuc9/sif/scmocha_latest.sif"
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
          "nohup java -Dconfig.file=/home/liuc9/github/scMOCHA/config/singularity.slurm.conf \\",
          "-jar /home/liuc9/tools/cromwell-78.jar \\",
          "run /home/liuc9/github/scMOCHA/scMOCHA-image.wdl \\",
          "-i {.jsonfile} 1>{.logfile} 2>{.errfile} &" |> glue::glue()
        )
        
        readr::write_lines(
          x = runwdl_cmd,
          file = runwdl_sh_file
        )
        
        
        
        
        .x |>
          dplyr::mutate(
            bname = gsub(
              pattern = "(.*)_L00(.*)_(.*)_001.fastq.gz",
              replacement = "\\3",
              x = name
            )
          ) |> 
          dplyr::mutate(
            targetname = glue::glue("{.y}_S1_L001_{bname}_001.fastq.gz")
          ) |> 
          dplyr::mutate(
            targetpath = file.path(.ydir, targetname)
          ) ->
          .xd
        # .xd$path
        # .xd$targetpath
        .xd |> 
          dplyr::mutate(
            a = purrr::map2(
              .x = path,
              .y = targetpath,
              .f = \(.p, .tp) {
                file.symlink(
                  from = .p,
                  to = .tp
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
  x = conf_scmocha |> dplyr::select(-data),
  file = file.path(
    "/home/liuc9/github/scMOCHA/runs/03-ADKP",
    "runfiles.csv"
  )
)

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
