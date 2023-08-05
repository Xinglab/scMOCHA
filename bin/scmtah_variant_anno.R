# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Mon Apr 10 16:56:33 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(httr)

# args --------------------------------------------------------------------
# load(file = "/scr1/users/liuc9/tmp/mito/flu2-a/cromwell-executions/SCMTAH/0138fcd0-c384-42c2-8704-6647767610d2/call-plot_scmtah/execution/scmtah.rda")

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------


# body --------------------------------------------------------------------


# Get variants ------------------------------------------------------------

cell_raw_cluster_forplot$forplot |> 
  dplyr::filter(depth > 3) |> 
  dplyr::mutate(variant = gsub(
    pattern = "[0-9]*",
    replacement = "",
    x = variant
  )) |> 
  tidyr::separate(
    col = variant,
    into = c("ref", "var")
  ) |> 
  dplyr::select(
    sample = barcode, 
    pos = pos, 
    ref = ref,
    var = var
  ) ->
  cell_variants

# 
readr::write_tsv(
  x = cell_variants,
  file = "cell_snvlist.tsv"
)
# library(httr)

tryCatch(
  {
    response <- POST(
      "https://mitomap.org/mitomaster/websrvc.cgi",
      body = list(
        file = upload_file("cell_snvlist.tsv"),
        fileType = "snvlist",
        output = "detail"
      ),
      encode = "multipart"
    )
    print(content(response, "text"))
  },
  error = function(err) {
    print(paste("HTTP error:", err$message))
  },
  warning = function(w) {
    print(paste("Warning:", w$message))
  },
  finally = {
    print("Done.")
  }
)


cluster_cluster_forplot$forplot |> 
  dplyr::filter(depth > 3) |>  
  dplyr::select(pos, variant) |>  
  dplyr::distinct() |> 
  dplyr::mutate(variant = gsub(
    pattern = "[0-9]*",
    replacement = "",
    x = variant
  )) |> 
  tidyr::separate(
    col = variant,
    into = c("ref", "var")
  ) |> 
  dplyr::mutate(sample = "Sample1") |> 
  dplyr::select(
    sample = sample, 
    pos = pos, 
    ref = ref,
    var = var
  ) ->
  cluster_variants

readr::write_tsv(
  x = cluster_variants,
  file = "cluster_snvlist.tsv"
)
tryCatch(
  {
    response <- POST(
      "https://mitomap.org/mitomaster/websrvc.cgi",
      body = list(
        file = upload_file("cluster_snvlist.tsv"),
        fileType = "snvlist",
        output = "detail"
      ),
      encode = "multipart"
    )
    # print(content(response, "text"))
  },
  error = function(err) {
    print(paste("HTTP error:", err$message))
  },
  warning = function(w) {
    print(paste("Warning:", w$message))
  },
  finally = {
    print("Done.")
  }
)

a <- content(
  x = response,
  as = "text",
  encoding = "UTF-8"
) |> 
  data.table::fread(
    sep = "\t"
  )

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------