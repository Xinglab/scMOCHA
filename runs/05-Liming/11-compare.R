datadir <- "/mnt/isilon/u01_project/PT/Comparison_from_different_cutoff_combinations_0708/Comparison_result_500_8000_2nd_C0.05_S0.2"

tibble::tibble(
  path = list.files(
    datadir,
    pattern = ".txt",
    full.names = T
  )
) |> 
  dplyr::mutate(
    d = purrr::map(
      .x = path,
      .f = data.table::fread
    )
  ) |> 
  dplyr::mutate(
    n = basename(path)
  ) ->
  d


d |> 
  dplyr::select(n, d) |> 
  dplyr::mutate(n = gsub(".txt", "", n)) |> 
  tidyr::separate(
    col = n,
    into = c("group", "mm", "celltype")
  ) |> 
  dplyr::select(group, celltype, d) |> 
  dplyr::mutate(nv = purrr::map_int(d, nrow)) |> 
  # dplyr::mutate(nv = c(5,10,22,5,50,21,60,23,8,69,14,72)) |> 
  dplyr::mutate(
    group2 = plyr::revalue(
      x = group,
      replace = c(
        "Chunjie" = "scRNASeq",
        "common" = "Shared",
        "Shiping" = "PacBio"
      )
    )
  ) |> 
  dplyr::mutate(
    group2 = factor(
      group2,
      levels = c("scRNASeq", "Shared", "PacBio")
    )
  ) ->
  dd

dd |> 
  dplyr::group_by(celltype) |> 
  dplyr::mutate(ratio = nv / sum(nv)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(
    label = glue::glue("{nv}\n({round(ratio, 3) * 100}%)")
  ) ->
  ddd

ddd |> 
  ggplot(aes(
    x = celltype,
    y = group2
  )) +
  geom_tile(
    aes(fill = ratio)
  ) +
  geom_text(
    aes(label = label),
    size = 8
  ) +
  scale_fill_gradient(
    low = "white",
    high = "red"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 18
    ),
    legend.position = "none ",
  ) ->
  p_share;p_share
ggsave(
  filename = "scRNASeq-PacBio-comparison.pdf",
  plot = p_share,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth2",
  height = 7,
  width = 10
)



dd |> 
  dplyr::filter(celltype == "143B") |> 
  dplyr::filter(group != "common") ->
  m
m$d[[1]]$Variant
m$d[[2]]$Variant


tibble::tibble(
  path = list.files(
    "/scr1/users/liuc9/mitochondrial/realdata/05-Liming/scmocha-mixed-cellline-high-depth/cromwell-executions/scMOCHA/8b0486c9-def1-40ed-896f-ea6dc9be39ca/call-gather_outputfiles/execution/WT",
    "*cluster.*.txt.gz*",
    full.names = T
  )
) |> 
  dplyr::mutate(d = purrr::map(path, data.table::fread)) |> 
  dplyr::mutate(n = basename(path)) |> 
  dplyr::mutate(n = gsub("cluster.|.txt.gz", "", n)) |> 
  dplyr::select(n, d) ->
  cluster



cluster |> 
  dplyr::filter(n != "coverage") |> 
  tidyr::unnest(cols = d) |> 
  dplyr::mutate(nv = V3 + V4) |> 
  dplyr::filter(V2 != "cluster_4") |> 
  dplyr::select(gt = n, pos = V1, group = V2, nv) ->
  cluster_n

cluster_n |> 
  dplyr::filter(pos == 2617) |> 
  dplyr::mutate(gt = factor(gt, levels = c("A", "G", "C", "T"))) |> 
  dplyr::group_by(gt) |> 
  dplyr::mutate(ratio = nv / sum(nv)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(
    label = glue::glue("{nv}\n({round(ratio, 3) * 100}%)")
  ) |> 
  dplyr::mutate(
    group2 = plyr::revalue(x = group, replace = c("cluster_0" = "A549", "cluster_1" = "WAL2A", "cluster_2" = "143B", "cluster_3" = "HEK293"))
  ) ->
  cluster_nn

cluster_nn |> 
  ggplot(aes(x = gt, y = group2)) +
  geom_tile(
    aes(fill = ratio)
  ) +
  geom_text(aes(label = label)) +
  scale_fill_gradient(
    low = "white",
    high = "red"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 18
    ),
    legend.position = "none ",
    plot.title = element_text(
      size = 16,
      hjust = 0.5
    )
  ) +
  labs(
    title = "9494A>G"
  )


mt_exons_df <- "/scr1/users/liuc9/mitochondrial/realdata/05-Liming/scmocha-mixed-cellline-high-depth/cromwell-executions/scMOCHA/023d7328-9097-4e50-8c11-19f860c5519e/call-cellranger_count/inputs/2014965526/mt_exons.df.rds.gz"
gtf_gene_df <-
  readr::read_rds(
    file = mt_exons_df
  )
library(gggenes)
ggplot(gtf_gene_df, aes(xmin = start, xmax = end, y = seqnames)) +
  # geom_gene_arrow() +
  geom_gene_arrow(
    aes(
      fill = gene_biotype
    ),
    arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")
  ) +
  scale_fill_brewer(
    palette = "Set1",
    name = "Gene type",
    labels = c("MT rRNA", "MT tRNA", "Protein coding")
  ) +
  ggrepel::geom_text_repel(
    aes(x = (start + end) / 2, label = gene_name, color = gene_biotype),
    # fill = "white",
    # nudge_x =1,
    # nudge_y = -0.1,
    size = 3,
    show.legend = F,
    max.overlaps = Inf,
  ) +
  scale_color_brewer(palette = "Set1") + 
  scale_x_continuous(
    limits = c(0, 17000),
    breaks = seq(0, 17000, 1000),
    expand = expansion(mult = c(0, 0.03)),
  ) +
  scale_y_discrete(
    expand = expansion(mult = c(0, 0), add = c(0, 0))
  ) +
  theme_genes() +
  theme(
    legend.position = "bottom",
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 14),
    legend.text = element_text(size = 14)
  ) ->
  pg;pg


cluster |> 
  dplyr::filter(n == "coverage") |> 
  tidyr::unnest(cols = d) |> 
  dplyr::select(pos = V1, group = V2, nv = V3) |> 
  dplyr::filter(group != "cluster_4") |> 
  dplyr::mutate(
    group2 = plyr::revalue(x = group, replace = c("cluster_0" = "A549", "cluster_1" = "WAL2A", "cluster_2" = "143B", "cluster_3" = "HEK293"))
  ) ->
  cluster_m

thepos <- c(302, 310, 4104, 4769, 7028, 7256, 7521, 8860)

cluster_m |> 
  ggplot(aes(x = pos, y = nv)) +
  geom_line(aes(color = group2 )) +
  geom_vline(xintercept = thepos, color = "red") +
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0)),
    limits = c(1, 17000),
    breaks = seq(0, 17000, 2000),
    labels = seq(0, 17000, 2000)
  ) +
  # ggsci::scale_color_jco(
  scale_color_brewer(
    name = "Dataset",
    palette = "Set1"
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    # axis.line.x.bottom = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = c(0.8, 0.2),
    legend.key = element_blank(),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    legend.text = element_text(
      size = 14,
      color = "black"
    ),
    legend.title = element_text(
      size = 16,
      colour = "black"
    )
  ) +
  labs(y = "Depth") ->
  pp

wrap_plots(
  pp,
  pg,
  ncol = 1,
  heights = c(0.9, 0.1)
) ->
  pg_merged_p_read_depth;pg_merged_p_read_depth

thepos <- c(302, 4104, 4769, 7028, 7256, 7521, 8860)
cluster_m |> 
  dplyr::filter(pos %in% thepos) |> 
  dplyr::mutate(pos = as.character(pos)) |> 
  ggplot(aes(x = pos, y = group2)) +
  geom_tile(aes(fill = nv)) +
  geom_text(aes(label = nv)) +
  scale_fill_gradient(
    low = "white",
    high = "red"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 18
    ),
    legend.position = "none ",
    plot.title = element_text(
      size = 16,
      hjust = 0.5
    )
  )
fasta <- Biostrings::readDNAStringSet("/home/liuc9/github/scMOCHA/fasta/rCRS.chrM.fasta")

fasta$chrM |> as.data.frame() |> 
  tibble::rownames_to_column(var = "pos") |> 
  dplyr::rename(ref = x) |> 
  dplyr::mutate(posref = glue::glue("{pos}{ref}")) |> 
  dplyr::mutate(pos = as.integer(pos)) ->
  fasta_df

cluster_n |> 
  dplyr::left_join(fasta_df, by = "pos") |> 
  # dplyr::mutate(pos = as.character(pos)) |> 
  dplyr::mutate(gt = factor(gt, levels = c("A", "G", "C", "T"))) |> 
  dplyr::group_by(group, pos) |> 
  # dplyr::group_by(pos, gt) |> 
  dplyr::mutate(ratio = nv / sum(nv)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(
    label = glue::glue("{nv}\n({round(ratio, 3) * 100}%)")
  ) |> 
  dplyr::mutate(
    group2 = plyr::revalue(x = group, replace = c("cluster_0" = "A549", "cluster_1" = "WAL2A", "cluster_2" = "143B", "cluster_3" = "HEK293"))
  ) ->
  cluster_n_forplot

cluster_n_forplot |> 
  dplyr::filter(pos %in% thepos) |> 
  dplyr::mutate(pos = as.character(pos)) |> 
  ggplot(aes(x = posref, y = gt)) +
  geom_tile(aes(fill = nv)) +
  geom_text(aes(label= label), size = 3.5) +
  scale_fill_gradient(
    low = "white",
    high = "red"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    # axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 18
    ),
    legend.position = "none ",
    plot.title = element_text(
      size = 16,
      hjust = 0.5
    ),
    strip.background = element_rect(
      fill = NA,
      color = "black",
    ),
    strip.text = element_text(
      color = "black",
      size = 14,
      face = "bold"
    ),
    axis.line = element_line(
      color = "black"
    )
  ) +
  facet_wrap(~group2, ncol = 1, strip.position = "right") ->
  p_tile;p_tile

ggsave(
  filename = "PacBio-specific-143B.pdf",
  plot = p_tile,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth2",
  height = 6,
  width = 9
)

dd |> 
  dplyr::filter(group != "common") |> 
  dplyr::filter(celltype == "143B") |> 
  tidyr::unnest(cols = d) |> 
  dplyr::select(v = Variant, group2) |> 
  dplyr::mutate(
    pos = gsub(
      ">|[AGCT]",
      "", x = v
    )
  ) |> 
  dplyr::select(pos, tech = group2) |> 
  dplyr::mutate(
    pos = as.integer(pos),
    tech = as.character(tech)
  ) |> 
  dplyr::distinct(pos) |> 
  dplyr::arrange(pos) ->
  c_143B

cluster_n_forplot |> 
  dplyr::inner_join(c_143B, by = "pos") |> 
  dplyr::arrange(pos) |> 
  dplyr::mutate(
    posref = forcats::fct_reorder(
      .f = posref,
      .x = pos
    )
  ) |> 
  ggplot(aes(x = posref, y = gt)) +
  geom_tile(aes(fill = nv)) +
  geom_text(aes(label= label), size = 3.5) +
  scale_fill_gradient(
    low = "white",
    high = "red"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    # axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 18
    ),
    legend.position = "none ",
    plot.title = element_text(
      size = 16,
      hjust = 0.5
    ),
    strip.background = element_rect(
      fill = NA,
      color = "black",
    ),
    strip.text = element_text(
      color = "black",
      size = 14,
      face = "bold"
    ),
    axis.line = element_line(
      color = "black"
    )
  ) +
  facet_wrap(~group2, ncol = 1, strip.position = "right") ->
  p_tile;p_tile

ggsave(
  filename = "PacBio-specific-143B.pdf",
  plot = p_tile,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth2",
  height = 6,
  width = 11
)



dd |> 
  dplyr::filter(group != "common") |> 
  dplyr::filter(celltype == "HEK293") |> 
  tidyr::unnest(cols = d) |> 
  dplyr::select(v = Variant, group2) |> 
  dplyr::mutate(
    pos = gsub(
      ">|[AGCT]",
      "", x = v
    )
  ) |> 
  dplyr::select(pos, tech = group2) |> 
  dplyr::mutate(
    pos = as.integer(pos),
    tech = as.character(tech)
  ) |> 
  dplyr::distinct(pos, .keep_all = T) |> 
  dplyr::arrange(pos) ->
  c_HEK293

cluster_n_forplot |> 
  dplyr::inner_join(
    c_HEK293 |> 
      dplyr::filter(tech == "scRNASeq"), 
    by = "pos"
  ) |> 
  dplyr::arrange(pos) |> 
  dplyr::mutate(
    posref = forcats::fct_reorder(
      .f = posref,
      .x = pos
    )
  ) |> 
  ggplot(aes(x = posref, y = gt)) +
  geom_tile(aes(fill = nv)) +
  geom_text(aes(label= label), size = 3.5) +
  scale_fill_gradient(
    low = "white",
    high = "red"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    # axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 18
    ),
    legend.position = "none ",
    plot.title = element_text(
      size = 16,
      hjust = 0.5
    ),
    strip.background = element_rect(
      fill = NA,
      color = "black",
    ),
    strip.text = element_text(
      color = "black",
      size = 14,
      face = "bold"
    ),
    axis.line = element_line(
      color = "black"
    )
  ) +
  facet_wrap(~group2, ncol = 1, strip.position = "right") ->
  p_tile;p_tile

ggsave(
  filename = "PacBio-specific-HEK293-scRNASeq.pdf",
  plot = p_tile,
  device = "pdf",
  path = "/home/liuc9/github/scMOCHA/05-Liming/scmocha-mixed-cellline-high-depth2",
  height = 6,
  width = 19
)
