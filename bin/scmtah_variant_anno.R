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
  cluster_variants




# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------