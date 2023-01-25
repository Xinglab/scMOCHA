# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Jan 24 15:18:16 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
# library(ggcoverage)
# 
# meta.file <- system.file("extdata", "RNA-seq", "meta_info.csv", package = "ggcoverage")
# sample.meta = read.csv(meta.file)
# 
# # track folder
# track.folder = system.file("extdata", "RNA-seq", package = "ggcoverage")
# # load bigwig file
# track.df = ggcoverage::LoadTrackFile(track.folder = track.folder, format = "bw",
#                          meta.info = sample.meta)
# body --------------------------------------------------------------------

coverage <- readr::read_tsv(
  file = "/scr1/users/liuc9/mitochondrial/PBMC_10k_v3_10x/MT.coverage",
  col_names = c("chr", "pos", "depth")
)



gtf107 <- readr::read_rds(file = "~/tmp/Homo_sapiens.GRCh38.107.gtf.plyranges.rds")


library(ggtranscript)

gtf_gene <- gtf107 %>% 
  plyranges::filter(seqnames == "MT") %>% 
  plyranges::filter(type == "exon")

gtf_gene %>% 
  as.data.frame() ->
  gtf_gene_df

coverage %>%
  ggplot(aes(x=pos, y = depth)) +
  # geom_line()
  geom_bar(stat = "identity") +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.03)),
    limits = c(1, 17000),
  ) +
  scale_y_continuous(
    expand = c(0.01, 0)
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  labs(y = "Depth") ->
  p1;p1

gtf_gene_df %>% 
  ggplot(aes(
    xstart = start,
    xend = end,
    y = gene_name
  )) +
  geom_range( aes(fill = transcript_biotype)) +
  geom_intron(
    data = to_intron(gtf_gene_df, "transcript_name"),
    aes(strand = strand)
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.03)),
    limits = c(1, 17000),
    breaks = seq(1000, 17000, 1000),
    labels = seq(1000, 17000, 1000)
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    axis.line = element_line(color = "black"),
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.title.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_blank(),
    legend.position = "bottom"
    ) +
  labs(
    x = "Position",
  ) ->
  p2;p2

p <- cowplot::plot_grid(
  plotlist = list(p1, p2),
  ncol = 1,
  align = "v",
  rel_heights = c(0.4, 0.6)
)

ggsave(
  filename = "reads-coverage.pdf",
  plot = p,
  device = "pdf",
  path = "data/PBMC_10k_v3_10x/result/03-coverage",
  width = 12,
  height = 9
)

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "data/PBMC_10k_v3_10x/rda/coverage.rda")
