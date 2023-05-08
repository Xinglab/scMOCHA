library(pasillaBamSubset)

library(gmoviz)

ideogram_1 <- GRanges(
  seqnames = c("chrA", "chrB", "chrC"),
  ranges = IRanges(start = rep(0, 3), end = rep(1000, 3))
)

ideogram_2 <- data.frame(
  chr = c("chrA", "chrB", "chrC"),
  start = rep(0, 3),
  end = rep(1000, 3)
)


fly_ideogram <- getIdeogramData(bam_file = pasillaBamSubset::untreated3_chr4())

fly_ideogram_chr4_only <- getIdeogramData(
  fasta_file = pasillaBamSubset::dm3_chr4()
)

mt_ideogram <- getIdeogramData(
  fasta_file = "/home/liuc9/github/scMOCHA/fasta/rCRS.MT.fasta"
)


gmovizInitialise(fly_ideogram_chr4_only, track_height = 0.15)


gmovizInitialise(
  mt_ideogram,
  track_height = 0.15
)

gmovizInitialise(
  fly_ideogram_chr4_only,
  space_between_sectors = 25, # bigger space between start & end
  start_degree = 78, # rotate the circle
  sector_label_size = 1, # bigger label
  track_height = 0.15, # thicker rectangle
  xaxis_spacing = 30
)

mt_ideogram <- getIdeogramData(
  fasta_file = "/home/liuc9/github/scMOCHA/fasta/rCRS.MT.fasta"
)
gmovizInitialise(
  mt_ideogram,
  space_between_sectors = 25, # bigger space between start & end
  start_degree = 78, # rotate the circle
  sector_label_size = 1, # bigger label
  track_height = 0.15, # thicker rectangle
  xaxis_spacing = 30
)


chr4_coverage <- getCoverage(
  regions_of_interest = "chr4",
  bam_file = pasillaBamSubset::untreated3_chr4(),
  window_size = 100
)

mt_mono <- getCoverage(
  regions_of_interest = "MT",
  bam_file = "/scr1/users/liuc9/tmp/mito/flu2-a/cromwell-executions/SCMTAH/8d90856b-1c1f-48a1-823c-b38ed32ba82d/call-cell_cluster_annotation/execution/MT_cluster.TAG_CJ_Mono.bam",
  window_size = 100
)

mt_mono |>
  data.table::as.data.table() |>
  dplyr::mutate(seqnames = "MT") |>
  plyranges::as_granges() ->
mt_mono_a

gmovizInitialise(
  ideogram_data = fly_ideogram_chr4_only,
  coverage_data = chr4_coverage,
  coverage_rectangle = "chr4",
  xaxis_spacing = 30
)

gmovizInitialise(
  ideogram_data = mt_ideogram,
  coverage_data = mt_mono_a,
  coverage_rectangle = "MT",
  xaxis_spacing = 30
)

label <- GRanges(
  seqnames = "chr4",
  ranges = IRanges(start = 240000, end = 280000),
  label = "region A"
)
gmovizInitialise(
  fly_ideogram_chr4_only,
  label_data = label,
  space_between_sectors = 25,
  start_degree = 78,
  sector_label_size = 1,
  xaxis_spacing = 30,
  coverage_data = chr4_coverage,
  coverage_rectangle = "chr4"
)

labels_from_file <- getLabels(
  gff_file = system.file("extdata", "example.gff3", package = "gmoviz"),
  colour_code = TRUE
)
gmovizInitialise(
  fly_ideogram_chr4_only,
  label_data = labels_from_file,
  label_colour = labels_from_file$colour,
  space_between_sectors = 25, start_degree = 78,
  sector_label_size = 1, xaxis_spacing = 30,
  coverage_data = chr4_coverage,
  coverage_rectangle = "chr4"
)


features <- getFeatures(
  gff_file = system.file("extdata", "example.gff3", package = "gmoviz"),
  colours = gmoviz::rich_colours
)
gmovizInitialise(
  fly_ideogram_chr4_only,
  space_between_sectors = 25,
  start_degree = 78, xaxis_spacing = 30,
  sector_label_size = 1
)
drawFeatureTrack(
  features,
  feature_label_cutoff = 80000,
  track_height = 0.18
)


numeric_data <- GRanges(
  seqnames = rep("chr4", 50),
  ranges = IRanges(
    start = sample(0:1320000, 50),
    width = 1
  ),
  value = runif(50, 0, 25)
)
gmovizInitialise(
  fly_ideogram_chr4_only,
  space_between_sectors = 25, start_degree = 78,
  sector_label_size = 1, xaxis_spacing = 30
)
## scatterplot
# drawScatterplotTrack(numeric_data)

## line graph
drawLinegraphTrack(
  sort(numeric_data),
  gridline_colour = NULL
)
# drawLinegraphTrack(
#   sort(numeric_data),
#   gridline_colour = NULL
# )
# drawLinegraphTrack(
#   sort(numeric_data),
#   gridline_colour = NULL
# )
drawFeatureTrack(features, feature_label_cutoff = 80000, track_height = 0.15)


# mito --------------------------------------------------------------------



mt_ideogram <- getIdeogramData(
  fasta_file = "/home/liuc9/github/scMOCHA/fasta/rCRS.MT.fasta"
)
gmovizInitialise(
  mt_ideogram,
  space_between_sectors = 25, # bigger space between start & end
  start_degree = 78, # rotate the circle
  sector_label_size = 1, # bigger label
  track_height = 0.15, # thicker rectangle
  xaxis_spacing = 30
)



mt_mono <- getCoverage(
  regions_of_interest = "MT",
  bam_file = "/scr1/users/liuc9/tmp/mito/flu2-a/cromwell-executions/SCMTAH/8d90856b-1c1f-48a1-823c-b38ed32ba82d/call-cell_cluster_annotation/execution/MT_cluster.TAG_CJ_Mono.bam",
  window_size = 50
)

library(Rsamtools)
bam <- BamFile("/scr1/users/liuc9/tmp/mito/flu2-a/cromwell-executions/SCMTAH/8d90856b-1c1f-48a1-823c-b38ed32ba82d/call-cell_cluster_annotation/execution/MT_cluster.TAG_CJ_other.bam")

# use the indexBam function to create an index file
indexBam(bam)

mt_mono <- getCoverage(
  regions_of_interest = "MT",
  bam_file = "/scr1/users/liuc9/tmp/mito/flu2-a/cromwell-executions/SCMTAH/8d90856b-1c1f-48a1-823c-b38ed32ba82d/call-cell_cluster_annotation/execution/MT_cluster.TAG_CJ_other.bam",
  window_size = 50
)

mt_mono |>
  data.table::as.data.table() |>
  dplyr::mutate(seqnames = "MT") |>
  plyranges::as_granges() ->
mt_mono_a

# conn <- DBI::dbConnect(
#   duckdb::duckdb(),
#   "/home/liuc9/data/refdata/ensembl/Homo_sapiens.GRCh38.107.gtf.plyranges.duckdb"
# )
# dplyr::tbl(conn, "grch38_107_plyranges") |>
#   dplyr::filter(seqnames == "MT") |>
#   dplyr::filter(type == "gene") |>
#   dplyr::mutate(
#     label = gene_name,
#     track = 1,
#     type = gene_biotype
#   ) |>
#   data.table::as.data.table() |>
#   dplyr::mutate(
#     colour = dplyr::case_match(
#       type,
#       "Mt_tRNA" ~ ggsci::pal_aaas()(4)[1],
#       "Mt_rRNA" ~ ggsci::pal_aaas()(4)[2],
#       "protein_coding" ~ ggsci::pal_aaas()(4)[3]
#     )
#   ) |>
#   dplyr::mutate(
#     shape = ifelse(
#       type == "Mt_tRNA",
#       "rectangle",
#       "forward_arrow"
#     )
#   ) |>
#   dplyr::select(
#     seqnames, start, end, width, strand,
#     label, track, type, shape, colour
#   ) |>
#   plyranges::as_granges() ->
#   mt_features
# DBI::dbDisconnect(conn, shutdown=TRUE)

# readr::write_rds(
#   mt_features,
#   "/home/liuc9/github/scMOCHA/fasta/mt_features.grange.gmoviz.rds.gz",
#   compress = "gz"
# )
mt_features <- readr::read_rds(
  "/home/liuc9/github/scMOCHA/fasta/mt_features.grange.gmoviz.rds.gz"
)


mt_features |>
  plyranges::filter(type == "Mt_tRNA") ->
  mt_features_pc
mt_features |>
  plyranges::filter(type != "Mt_tRNA") ->
  mt_features_npc

gmovizInitialise(
  mt_ideogram,
  space_between_sectors = 25,
  start_degree = 78,
  xaxis_spacing = 30,
  sector_label_size = 1,
  coverage_data = mt_mono_a,
  coverage_rectangle = "MT"
)
drawFeatureTrack(
  mt_features_npc,
  # feature_label_cutoff = 80000,
  track_height = 0.13
)
drawFeatureTrack(
  mt_features_pc,
  feature_label_cutoff = 80000,
  # track_height = 0.18
)



