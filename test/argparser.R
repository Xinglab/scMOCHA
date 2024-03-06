#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Mar  5 13:37:28 2024
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
library(data.table)
#library(rlang)

# args --------------------------------------------------------------------

'Naval Fate.

Usage:
  naval_fate.R ship new <name>...
  naval_fate.R ship <name> move <x> <y> [--speed=<kn>]
  naval_fate.R ship shoot <x> <y>
  naval_fate.R mine (set|remove) <x> <y> [--moored | --drifting]
  naval_fate.R (-h | --help)
  naval_fate.R --version

Options:
  -h --help     Show this screen.
  --version     Show version.
  --speed=<kn>  Speed in knots [default: 10].
  --moored      Moored (anchored) mine.
  --drifting    Drifting mine.

' -> doc

library(docopt)
arguments <- docopt(doc, version = 'Naval Fate 2.0')
print(arguments)

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------


# body --------------------------------------------------------------------



# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------