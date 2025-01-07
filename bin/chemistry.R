#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Sep 12 11:58:54 2024
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
library(prismatic)
library(paletteer)
library(data.table)
# library(rlang)
library(GetoptLong)
library(logger)

# args --------------------------------------------------------------------
# chemistry_detector <- "/home/liuc9/github/scMOCHA-data/data/GSE181279/cromwell-executions/scMOCHABatch/2a7602ab-7f25-4b55-8595-58fa7039c5d8/call-scMOCHA/shard-0/sub.scMOCHA/65c88c67-7a86-4c32-935d-fd9361e7a5b1/call-cellranger_count/execution/GSM5494107/SC_RNA_COUNTER_CS/SC_MULTI_CORE/MULTI_CHEMISTRY_DETECTOR/DETECT_COUNT_CHEMISTRY/fork0/_outs"
# output_id <- "GSM5494107"

# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.
verbose <- FALSE
spec <- "
Usage: Rscript foorbar.R [options]

Options:
<output_id=s> output id
<verbose!> Print messages
"

GetoptLong.options(help_style = "two-column")
GetoptLong(spec, template_control = list(opt_width = 21))

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------
log_threshold(TRACE)
log_layout(layout_glue_colors)

# log_info('Starting the script...')
# log_debug('This is the second log line')
# log_trace('Note that the 2nd line is being placed right after the 1st one.')
# log_success('Doing pretty well so far!')
# log_warn('But beware, as some errors might come :/')
# log_error('This is a problem')
# log_debug('Note that getting an error is usually bad')
# log_error('This is another problem')
# log_fatal('The last problem')

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
chemistry_detector <- file.path(
  output_id,
  "SC_RNA_COUNTER_CS/SC_MULTI_CORE/MULTI_CHEMISTRY_DETECTOR/DETECT_COUNT_CHEMISTRY/fork0/_outs"
)
chem <- jsonlite::fromJSON(chemistry_detector)

# body --------------------------------------------------------------------

chem_description <- chem$chemistry_defs$`Gene Expression`$description

chem_endedness <- chem$chemistry_defs$`Gene Expression`$endedness

chem_name <- chem$chemistry_defs$`Gene Expression`$name

chem_df <- tibble::tibble(
  name = chem_name,
  endedness = chem_endedness,
  description = chem_description
)

# output ------------------------------------------------------------------


outfile <- file.path(
  output_id,
  "outs",
  "chemistry.csv"
)

readr::write_csv(
  x = chem_df,
  file = outfile
)

readr::write_lines(
  x = chem_name,
  path = file.path(
    output_id,
    "outs",
    "chemistry_name.txt"
  )
)

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
