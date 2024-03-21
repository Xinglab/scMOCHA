#!/usr/bin/env Rscript
# pkgs <- c()


# pkgs <- unique(pkgs)
# ap.db <- available.packages(contrib.url(BiocManager::repositories()))
# ap <- rownames(ap.db)
# installed_packages <- rownames(installed.packages())
# pkgs_to_install <- setdiff(pkgs[pkgs %in% ap], installed_packages)
# pkgs_to_install

# BiocManager::install(pkgs_to_install, update=FALSE, ask=FALSE)
# From github

# options(buildtools.check = function(action) TRUE )
# devtools::install_github('immunogenomics/presto', ref = 'master', upgrade = "never", force = TRUE)
# devtools::install_github('satijalab/azimuth', ref = 'master', upgrade = "never", force = TRUE)
# devtools::install_github('satijalab/seurat-data', ref = 'master', upgrade = "never", force = TRUE)
# devtools::install_github("dzhang32/ggtranscript", ref = 'master', upgrade = "never", force = TRUE)

# suppressWarnings(BiocManager::install(update=TRUE, ask=FALSE))
# Remove tmp directory
#  SeuratData

# sd <- SeuratData::AvailableData() |>
#   dplyr::filter(grepl("Azimuth Reference", x = Summary)) |>
#   dplyr::filter(species == "human") |>
#   dplyr::filter(!Installed)

# purrr::map(
#   sd$Dataset,
#   \(.x) {
#       tryCatch(
#         expr = {
#           SeuratData::InstallData(.x)
#         },
#         error = function(e) {
#           1
#         }
#       )
#   }
# )


# just in case there were warnings, we want to see them
# without having to scroll up:
# warnings()


# if (!is.null(warnings()))
# {
#   w <- capture.output(warnings())
#   if (length(grep("is not available|had non-zero exit status", w))) quit(save="no", status=0L, runLast = FALSE)
# }


# unlink(x = '/tmp/*', recursive=T)