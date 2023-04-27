# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Apr 27 14:47:37 2023
# @DESCRIPTION: 


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
datadir <- "/scr1/users/liuc9/mitochondrial/realdata/01-Sci_Immunol_32651212"

sratable <- readr::read_csv(
  file = file.path(datadir, "SraRunTable.txt")
)

# body --------------------------------------------------------------------

sratable |> 
  dplyr::select(srrid = Run) |> 
  dplyr::mutate(
    srrdir = file.path(datadir, srrid)
  ) |> 
  dplyr::mutate(
    srrdir_exist = file.exists(srrdir)
  ) |> 
  dplyr::mutate(
    srafile = file.path(srrdir, glue::glue("{srrid}.sra"))
  ) |> 
  dplyr::mutate(
    srafile_exist = file.exists(srafile)
  ) ->
  srafiles



srafiles |> 
  dplyr::mutate(
    dump_slrm = purrr::map2_chr(
      .x = srrdir,
      .y = srafile,
      .f = function(.x, .y) {
        .srrid <- basename(.x)
        
        cmd_slrm <- c(
          "#!/usr/bin/env bash",
          "# @AUTHOR: Chun-Jie Liu",
          "# @CONTACT: chunjie.sam.liu.at.gmail.com",
          "# @DATE: {lubridate::now()}" |> glue::glue(),
          "",
          "#SBATCH --signal=USR2",
          "#SBATCH --ntasks=1",
          "#SBATCH --cpus-per-task=10",
          "#SBATCH --mem=50G",
          "#SBATCH --time=720:00:00",
          "#SBATCH --output={.x}/dump.job.{.srrid}.%j" |> glue::glue(),
          "module load R/4.1.0"
        )
        cmd_dump <- c(
          "/home/liuc9/tools/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump {.y} --include-technical --mem 50G --threads 10 --split-files --outdir {.x}" |> glue::glue()
        )
        cmd_rename <- c(
          "mv {.x}/{.srrid}_1.fastq {.x}/{.srrid}_S1_L001_R1_001.fastq" |> glue::glue(),
          "mv {.x}/{.srrid}_2.fastq {.x}/{.srrid}_S1_L001_R2_001.fastq" |> glue::glue()
        )
        
        cmd <- c(
          cmd_slrm,
          "",
          cmd_dump,
          "",
          cmd_rename
        )
        
        dump_slrm_file <-  file.path(
          .x,
          "dump_{.srrid}.slrm" |> glue::glue()
        )
        
        readr::write_lines(
          cmd,
          file = dump_slrm_file
        )
        dump_slrm_file
      }
    )
  ) ->
  sarfiles_dump


sarfiles_dump |> 
  dplyr::mutate(
    conf_json = purrr::map_chr(
      .x = srrdir,
      .f = function(.x) {
        .srrid <- basename(.x)
        
        conf <- list(
          "SCMTAH.output_id" = "{.srrid}" |> glue::glue(),
          "SCMTAH.fastqs" = "{.x}" |> glue::glue(),
          "SCMTAH.sample_id" = "{.srrid}" |> glue::glue(),
          "SCMTAH.transcriptome" = "/home/liuc9/data/refdata/mgatk_index/Human",
          "SCMTAH.rCRS" = "/home/liuc9/github/scRNAseq-MitoVariant/fasta/rCRS.MT.fasta",
          "SCMTAH.output_dir" = "{.srrid}" |> glue::glue()
        )
        
        conf_file <- file.path(
          .x,
          "{.srrid}.json" |> glue::glue()
        )
        
        jsonlite::write_json(
          x = conf,
          path = conf_file,
          auto_unbox = TRUE
        )
        
        conf_file
      }
    )
  ) ->
  sarfiles_dump_conf


sarfiles_dump_conf |> 
  dplyr::mutate(
    scmtah_sh = purrr::map_chr(
      .x = srrdir,
      .f = function(.x) {
        .srrid <- basename(.x)
        .jsonfile <- file.path(
          .x, "{.srrid}.json" |> glue::glue()
        )
        .errfile <- file.path(
          .x, "{.srrid}.err"
        )
        .logfile <- file.path(
          .x, "{.srrid}.log"
        )
        
        runwdl_sh_file <- file.path(
          .x, "runwdl_{.srrid}.sh" |> glue::glue()
        )
        
        runwdl_cmd <- c(
          "#!/usr/bin/env bash",
          "# @AUTHOR: Chun-Jie Liu",
          "# @CONTACT: chunjie.sam.liu.at.gmail.com",
          "# @DATE: {lubridate::now()}" |> glue::glue(),
          "",
          "module load Java/15.0.1",
          "nohup java -Dconfig.file=/home/liuc9/github/scRNAseq-MitoVariant/config/ref.conf \\",
          "-jar /home/liuc9/tools/cromwell-78.jar \\",
          "run /home/liuc9/github/scRNAseq-MitoVariant/scmtah.wdl \\",
          "-i {.jsonfile} 1>{.logfile} 2>{.errfile} &" |> glue::glue()
        )
        
        readr::write_lines(
          x = runwdl_cmd,
          file = runwdl_sh_file
        )
        
        runwdl_sh_file
      }
    )
  ) ->
  sarfiles_dump_conf_scmtah


readr::write_tsv(
  x = sarfiles_dump_conf_scmtah,
  file = "runs/01-Sci_Immunol_32651212/runfiles.tsv"
)

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------