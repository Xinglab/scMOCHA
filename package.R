pkgs <- c(
  ""
)


pkgs <- unique(pkgs)
ap.db <- available.packages(contrib.url(BiocManager::repositories()))
ap <- rownames(ap.db)
pkgs_to_install <- pkgs[pkgs %in% ap]

BiocManager::install(pkgs_to_install, update=FALSE, ask=FALSE)
# From github

remotes::install_github('satijalab/azimuth', ref = 'master', upgrade = "never", force = TRUE)
devtools::install_github('satijalab/seurat-data')

suppressWarnings(BiocManager::install(update=TRUE, ask=FALSE))
# Remove tmp directory

# just in case there were warnings, we want to see them
# without having to scroll up:
warnings()

if (!is.null(warnings()))
{
  w <- capture.output(warnings())
  if (length(grep("is not available|had non-zero exit status", w))) quit(save="no", status=0L, runLast = FALSE)
}

unlink(x = '/tmp/*', recursive=T)