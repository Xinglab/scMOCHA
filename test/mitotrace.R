library(MitoTrace)


bams <- list.files("/home/liuc9/tmp/bam_file", full.names = T, pattern = ".bam$")

fasta_loc <- "/scr1/users/liuc9/tmp/bam_file/GRCH38_MT.fa"


ann <- read.table("/scr1/users/liuc9/tmp/bam_file/SRR_donors_list")

mae_res <- MitoTrace(bam_list = bams, fasta = fasta_loc, chr_name = "MT")

mae_res


af <- calc_allele_frequency(mae_res)
colnames(af) <- gsub(".filter.bam", "", colnames(af))

af01 <- af[ ,na.omit(ann$V1[match(colnames(af), ann$V1)])]

colnames(af01) <- na.omit(ann$V2[match(colnames(af), ann$V1)])
donors <- unlist(lapply(colnames(af01), function(x) paste(strsplit(x, "_", fixed = T)[[1]][1:2], collapse = "_")))


colnames(af) <- unlist(lapply(colnames(af), function(x) strsplit(x, ".", fixed=T)[[1]][1]))


# plot first row
par(mfrow = c(2,6))
boxplot(rev(split(af01["779T>C", ], donors == "Donor1_C101")), main="779 T>C",ylab="Variant heteroplasmy", names=c("C101", "Other"))
boxplot(rev(split(af01["8978T>C", ], donors == "Donor1_C103")), main="8979 T>C", names=c("C103","Other"))
boxplot(rev(split(af01["6712A>G", ], donors == "Donor1_C107")), main="6712 A>G", names=c("C107","Other"))
boxplot(rev(split(af01["1082A>G", ], donors == "Donor1_C109")), main="1082 A>G", names=c("C109","Other"))
boxplot(rev(split(af01["3776G>A", ], donors == "Donor1_C112")), main="3776 T>A", names=c("C112","Other"))
boxplot(rev(split(af01["7275T>C", ], donors == "Donor1_C114")), main="7275 T>C", names=c("C114","Other"))
# plot second row
boxplot(rev(split(af01["13093G>A", ], donors == "Donor1_C116")), main="13093 G>A", ylab="Variant heteroplasmy", names=c("C116","Other"))
boxplot(rev(split(af01["7340G>A", ], donors == "Donor1_C118")), main="7340 G>A", names=c("C118","Other"))
boxplot(rev(split(af01["7754G>A", ], donors == "Donor1_C120")), main="7755 G>A", names=c("C120","Other"))
boxplot(rev(split(af01["2646G>A", ], donors == "Donor1_C124")), main="2648 G>A", names=c("C124","Other"))
boxplot(rev(split(af01["11622A>C", ], donors == "Donor1_C132")), main="11623 T>C", names=c("C132","Other"))
boxplot(rev(split(af01["1446A>G", ], donors == "Donor1_C135")), main="1448 A>G", names=c("C135","Other"))



mitotrace_parallel <- function(
    bam_list = bams, 
    fasta = fasta_loc, 
    chr_name = "MT",
    tag_name = "CB", 
    barcodes = NULL,
    min_read = 1000, 
    max_depth = 1e+06,
    min_base_quality = 25, 
    min_mapq = 30, 
    min_nucleotide_depth = 0,
    min_minor_allele_depth = 0
    ) {
  require(Rsamtools)
  require(Matrix)
  require(seqinr)
  bases <- c("A", "C", "G", "T")
  reffasta <- seqinr::read.fasta(fasta)
  mitoChr <- attr(reffasta, "name")
  maxpos <- length(reffasta[[1]])
  reference_genome <- data.frame(
    postion = 1:maxpos, 
    base = toupper(unname(reffasta)[[1]])[1:maxpos]
    )
  which <- GRanges(
    seqnames = chr_name, 
    ranges = IRanges(
    1,
    maxpos
  ))
  if (length(bam_list) == 1) {
    combinedBam <- TRUE
  }
  if (length(bam_list) > 1) {
    combinedBam <- FALSE
  }
  bam_name_list_array <- BamFileList(bam_list)
  if (combinedBam) {
    print("Extracting barcodes from single BAM file and running pileup for each barcode separately")
    params <- ScanBamParam(tag = tag_name, which = which)
    if (is.null(barcodes)) {
      barcodes <- scanBam(bam_list, param = params)
      good_barcodes <- names(which(table(barcodes[[1]][[1]][[1]]) >
        min_read))
    }
    if (!is.null(barcodes)) {
      good_barcodes <- barcodes
    }
    total_mpileups <- lapply(good_barcodes, function(x) {
      filter <- list(x)
      names(filter) <- tag_name
      pileup_bam <- pileup(bam_list, scanBamParam = ScanBamParam(
        tagFilter = filter,
        which = which
      ), pileupParam = PileupParam(
        distinguish_strands = FALSE,
        max_depth = max_depth, min_base_quality = min_base_quality,
        min_mapq = min_mapq, min_nucleotide_depth = min_nucleotide_depth,
        min_minor_allele_depth = min_minor_allele_depth
      ))
      bases <- c("A", "T", "C", "G")
      base_counts <- lapply(bases, function(base) {
        mutation <- subset(pileup_bam, nucleotide ==
          base)
        data.frame(mutation$pos, mutation$count)
      })
      names(base_counts) <- bases
      base_counts
    })
    names(total_mpileups) <- good_barcodes
  } else {
    print("Running pileup on each BAM file separately")
    bam_name_list_array <- BamFileList(bam_list)
    total_mpileups <- parallel::mclapply(
      bam_name_list_array, 
      function(x) {
      pileup_bam <- pileup(
        x,
        scanBamParam = ScanBamParam(which = which),
        pileupParam = PileupParam(
          distinguish_strands = FALSE,
          max_depth = max_depth, min_base_quality = min_base_quality,
          min_mapq = min_mapq, min_nucleotide_depth = min_nucleotide_depth,
          min_minor_allele_depth = min_minor_allele_depth
        )
      )
      
      bases <- c("A", "T", "C", "G")
      base_counts <- lapply(
        bases, 
        function(base) {
        mutation <- subset(pileup_bam, nucleotide ==
          base)
        data.frame(mutation$pos, mutation$count)
      })
      names(base_counts) <- bases
      base_counts
    },
    mc.cores = 10
    )
    nom <- unlist(lapply(bam_list, basename))
    names(total_mpileups) <- nom
  }
  
  
  
  res_counts <- lapply(
    bases, 
    function(base) {
    allpos <- 1:maxpos
    counts_base_allcells <- do.call(
      cbind, 
      parallel::mclapply(
      total_mpileups,
      function(x) {
        count_base_cell <- x[[base]]
        count_base_cell$mutation.count[match(
          allpos,
          count_base_cell$mutation.pos
        )]
      },
      mc.cores = 10
    ))
    counts_base_allcells[which(is.na(counts_base_allcells))] <- 0
    
    rownames(counts_base_allcells) <- paste0(
      allpos, reference_genome$base,
      ">", base
    )
    
    colnames(counts_base_allcells) <- names(total_mpileups)
    as(counts_base_allcells, "sparseMatrix")
  })
  
  names(res_counts) <- bases
  coverage <- matrix(0, nrow(res_counts[[1]]), ncol(res_counts[[1]]))
  
  rownames(coverage) <- as.character(1:maxpos)
  colnames(coverage) <- colnames(res_counts[[1]])
  
  
  lapply(1:4, function(x) coverage <<- coverage + data.matrix(res_counts[[x]]))
  
  res_counts2 <- lapply(bases, function(base) {
    tmp <- res_counts[[base]]
    tmp[which(reference_genome$base != base), ]
  })
  
  names(res_counts2) <- bases
  
  allpos <- unique(unlist(lapply(res_counts2, rownames)))
  pos_tmp <- allpos
  pos_tmp <- unlist(lapply(pos_tmp, function(x) {
    substr(
      x, 1,
      nchar(x) - 3
    )
  }))
  pos_tmp <- as.numeric(pos_tmp)
  allpos <- allpos[order(pos_tmp)]
  matr <- matrix(0, length(allpos), ncol(res_counts2[["A"]]))
  rownames(matr) <- allpos
  colnames(matr) <- colnames(res_counts2[["A"]])
  lapply(res_counts2, function(x) {
    ok <- intersect(rownames(matr), rownames(x))
    matr[ok, ] <<- data.matrix(x[ok, ])
  })
  counts <- matr
  counts <- as(counts, "sparseMatrix")
  coverage <- as(coverage, "sparseMatrix")
  return(list(read_counts = counts, coverage = coverage))
}


function (object) {
  tmp <- rownames(object[[1]])
  pos <- unlist(lapply(tmp, function(x) substr(
    x, 
    1,
    nchar(x) - 3)))
  pos <- as.numeric(pos)
  (object[[1]])/(object[[2]][pos, ] + 1e-04)
}