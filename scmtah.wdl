
workflow SCMTAH {
  # version of this pipeline
  String version = "CellRanger v7.0.1"

  String output_id = "cellranger"
  String fastqs
  String sample_id
  String transcriptome = "/home/liuc9/data/refdata/mgatk_index/Human"
  String output_dir

  String chrM = "MT"
  File rCRS = "/home/liuc9/github/scRNAseq-MitoVariant/fasta/rCRS.MT.fasta"

  String cellrefname = "pbmcref"
  String celllevel = "celltype.l1"


  # Runtime attributes
  String memory = "50 GB"
  Int boot_disk_size_gb = 12
  String disk_space = "50"
  Int cpu = 20

  parameter_meta {
      output_id: "Output ID"
      fastqs: "Directory of fastq files"
      sample_id: "Prefix name for fastq"
      transcriptome: "CellRanger-compatible transcriptome reference (can be generated with cellranger mkref)"
      expect_cells: "Expected number of recovered cells (defaults to 3000)"
      memory: "The minimum amount of RAM to use for the Cromwell VM"
      boot_disk_size_gb: "Size of disk (GB) where the docker image is booted by the Cromwell VM"
      disk_space: "Amount of disk space (GB) to allocate to the Cromwell VM"
      cpu: "The minimum number of cores to use for the Cromwell VM"
  }

  call cellranger_count {
      input:
        output_id = output_id,
        fastqs = fastqs,
        sample_id = sample_id,
        transcriptome = transcriptome,
        chrM = chrM,
        memory = memory,
        boot_disk_size_gb = boot_disk_size_gb,
        disk_space = disk_space,
        cpu = cpu
  }

  call cell_cluster_annotation {
    input:
      h5file = cellranger_count.filtered_feature_bc_matrix,
      mt_bam = cellranger_count.mt_bam,
      mt_bam_index = cellranger_count.mt_bam_index,
      refname = cellrefname,
      celllevel = celllevel,
      memory = memory,
      boot_disk_size_gb = boot_disk_size_gb,
      disk_space = disk_space,
      cpu = cpu
  }

  call call_mt_variants {
    input:
      sorted_bam = cellranger_count.sorted_bam,
      sorted_bam_index = cellranger_count.sorted_bam_index,
      barcodes_tsv = cellranger_count.barcodes_tsv,
      mt_bam = cellranger_count.mt_bam,
      mt_bam_index = cellranger_count.mt_bam_index,
      mt_cluster_bam = cell_cluster_annotation.mt_cluster_bam,
      mt_cluster_bam_index = cell_cluster_annotation.mt_cluster_bam_index,
      mt_bulk_bam = cell_cluster_annotation.mt_bulk_bam,
      mt_bulk_bam_index = cell_cluster_annotation.mt_bulk_bam_index,
      chrM = chrM,
      rCRS = rCRS,
      memory = memory,
      boot_disk_size_gb = boot_disk_size_gb,
      disk_space = disk_space,
      cpu = cpu
  }

  output {
      # version of this pipeline
      String pipeline_version = version

      # cellranger count
      File filtered_feature_bc_matrix = cellranger_count.filtered_feature_bc_matrix
      File metrics_summary = cellranger_count.metrics_summary
      File sorted_bam = cellranger_count.sorted_bam
      File sorted_bam_index = cellranger_count.sorted_bam_index
      File web_summary = cellranger_count.web_summary
      File barcodes = cellranger_count.barcodes
      File barcodes_tsv = cellranger_count.barcodes_tsv
      File features = cellranger_count.features
      File matrix = cellranger_count.matrix
      File mt_depth = cellranger_count.mt_depth
      File mt_depth_plot = cellranger_count.mt_depth_plot
      File mt_bam = cellranger_count.mt_bam
      File mt_bam_index = cellranger_count.mt_bam_index


      # cell_cluster_annotation
      File azimuth_rda = cell_cluster_annotation.azimuth_rda
      File barcode_cluster = cell_cluster_annotation.barcode_cluster
      File barcode_bulk = cell_cluster_annotation.barcode_bulk
      File celltype_ratio = cell_cluster_annotation.celltype_ratio
      File plot_metrics = cell_cluster_annotation.plot_metrics
      File plot_pie_celltype = cell_cluster_annotation.plot_pie_celltype
      File plot_qc = cell_cluster_annotation.plot_qc
      File plot_umap = cell_cluster_annotation.plot_umap
      File qc_cell_stats = cell_cluster_annotation.qc_cell_stats
      File sc_azimuth_rds_gz = cell_cluster_annotation.sc_azimuth_rds_gz
  }

  meta {
    author: "Chun-Jie Liu"
    email : "chunjie.sam.liu@gmail.com"
    description: "Analyze 3' single-cell RNA-seq data using the 10X v3 Genomics Cellranger pipeline."
  }
}

task cellranger_count {
    String output_id
    String fastqs
    String sample_id
    String transcriptome

    String chrM = "MT"

    String memory
    Int boot_disk_size_gb
    String disk_space
    Int cpu

    command {

      module load R/4.1.0

      # cell ranger output to get bams and h5 files
      cellranger count \
        --id=${output_id} \
        --fastqs=${fastqs} \
        --sample=${sample_id} \
        --transcriptome=${transcriptome} \
        --nosecondary \
        --disable-ui \
        --localcores ${cpu}

      # decompress barcodes.tsv.gz
      gunzip -c ${output_id}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > ${output_id}/outs/filtered_feature_bc_matrix/barcodes.tsv

      # extract MT.bam
      samtools view -b -h -o ${output_id}/outs/possorted_genome_bam.MT.bam ${output_id}/outs/possorted_genome_bam.bam ${chrM}
      samtools index ${output_id}/outs/possorted_genome_bam.MT.bam

      # mt depth
      samtools depth -a -r ${chrM} --threads=${cpu} ${output_id}/outs/possorted_genome_bam.MT.bam > ${output_id}/outs/possorted_genome_bam.MT.depth

      # Depth plot
      Rscript /home/liuc9/github/scRNAseq-MitoVariant/bin/depth.R ${output_id}/outs/possorted_genome_bam.MT.depth ${output_id}/outs/possorted_genome_bam.MT.depth.pdf

    }

    output {
        File filtered_feature_bc_matrix = "${output_id}/outs/filtered_feature_bc_matrix.h5"
        File metrics_summary = "${output_id}/outs/metrics_summary.csv"
        File sorted_bam = "${output_id}/outs/possorted_genome_bam.bam"
        File sorted_bam_index = "${output_id}/outs/possorted_genome_bam.bam.bai"
        File web_summary = "${output_id}/outs/web_summary.html"
        File barcodes = "${output_id}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
        File barcodes_tsv = "${output_id}/outs/filtered_feature_bc_matrix/barcodes.tsv"
        File features = "${output_id}/outs/filtered_feature_bc_matrix/features.tsv.gz"
        File matrix = "${output_id}/outs/filtered_feature_bc_matrix/matrix.mtx.gz"
        File mt_depth = "${output_id}/outs/possorted_genome_bam.MT.depth"
        File mt_depth_plot = "${output_id}/outs/possorted_genome_bam.MT.depth.pdf"
        File mt_bam = "${output_id}/outs/possorted_genome_bam.MT.bam"
        File mt_bam_index = "${output_id}/outs/possorted_genome_bam.MT.bam.bai"
    }
}

task cell_cluster_annotation {
  File h5file
  File mt_bam
  File mt_bam_index

  String refname = "pbmcref"
  String celllevel = "celltype.l1"

  String memory
  Int boot_disk_size_gb
  String disk_space
  Int cpu

  command {
    module load R/4.1.0
    # cell cluster annotation
    Rscript /home/liuc9/github/scRNAseq-MitoVariant/bin/azimuth.R ${h5file} ${refname} ${celllevel}

    # addtags for cluster
    sinto addtags \
      -b ${mt_bam} \
      -f barcode_cluster.tsv \
      -o MT_cluster.bam \
      -p ${cpu}
    samtools index MT_cluster.bam

    # addtags for bulk
    sinto addtags \
      -b ${mt_bam} \
      -f barcode_bulk.tsv \
      -o MT_bulk.bam \
      -p ${cpu}
    samtools index MT_bulk.bam

  }
  output {
    File azimuth_rda = "azimuth.rda"
    File barcode_cluster = "barcode_cluster.tsv"
    File barcode_bulk = "barcode_bulk.tsv"
    File celltype_ratio = "celltype_ratio.tsv"
    File plot_metrics = "plot-metrics.pdf"
    File plot_pie_celltype = "plot-pie-celltype.pdf"
    File plot_qc = "plot-qc.pdf"
    File plot_umap = "plot-umap.pdf"
    File qc_cell_stats = "qc-cell-stats.xlsx"
    File sc_azimuth_rds_gz = "sc_azimuth.rds.gz"
    File mt_cluster_bam = "MT_cluster.bam"
    File mt_cluster_bam_index = "MT_cluster.bam.bai"
    File mt_bulk_bam = "MT_bulk.bam"
    File mt_bulk_bam_index = "MT_bulk.bam.bai"
  }
}

task call_mt_variants {
  File sorted_bam
  File sorted_bam_index
  File barcodes_tsv

  File mt_bam
  File mt_bam_index
  File mt_cluster_bam
  File mt_cluster_bam_index
  File mt_bulk_bam
  File mt_bulk_bam_index

  String chrM
  File rCRS


  String memory
  Int boot_disk_size_gb
  String disk_space
  Int cpu

  command {

    module load R/4.1.0
    # call variants on single cell level
    mgatk tenx -i ${sorted_bam} \
      -o cell \
      -n cell \
      -g ${rCRS} \
      -c ${cpu} \
      -bt CB \
      -b ${barcodes_tsv} \
      -ub UB

    # the cell/final/ last "/" is important
    python /home/liuc9/github/scRNAseq-MitoVariant/bin/variant_calling_cell_raw.py \
      cell/final/ \
      cell \
      16569 \
      10 \
      ${chrM}

    # call variants on cluster level
    mgatk bcall -i ${mt_cluster_bam} \
      -o cluster \
      -n cluster \
      -g ${rCRS} \
      -c ${cpu} \
      -bt CJ

    python /home/liuc9/github/scRNAseq-MitoVariant/bin/variant_calling_cluster.py \
      cluster/final/ \
      cluster \
      16569 \
      10 \
      ${chrM}

  }
  output {
    # cell level
    File cell_A_txt_gz = "cell/final/cell.A.txt.gz"
    File cell_C_txt_gz = "cell/final/cell.C.txt.gz"
    File cell_G_txt_gz = "cell/final/cell.G.txt.gz"
    File cell_T_txt_gz = "cell/final/cell.T.txt.gz"
    File cell_cell_heteroplasmic_df_tsv_gz = "cell/final/cell.cell_heteroplasmic_df.tsv.gz"
    File cell_cell_heteroplasmic_df_raw_tsv_gz = "cell/final/cell.cell_heteroplasmic_df_raw.tsv.gz"
    File cell_coverage_txt_gz = "cell/final/cell.coverage.txt.gz"
    File cell_depthTable_txt = "cell/final/cell.depthTable.txt"
    File cell_rds = "cell/final/cell.rds"
    File cell_signac_rds = "cell/final/cell.signac.rds"
    File cell_variant_stats_tsv_gz = "cell/final/cell.variant_stats.tsv.gz"
    File cell_vmr_strand_plot_png = "cell/final/cell.vmr_strand_plot.png"

    # cluster level
    File barcodeQuants_tsv = "cluster/final/barcodeQuants.tsv"
    File cluster_A_txt_gz = "cluster/final/cluster.A.txt.gz"
    File cluster_C_txt_gz = "cluster/final/cluster.C.txt.gz"
    File cluster_G_txt_gz = "cluster/final/cluster.G.txt.gz"
    File cluster_T_txt_gz = "cluster/final/cluster.T.txt.gz"
    File cluster_cell_heteroplasmic_df_tsv_gz = "cluster/final/cluster.cell_heteroplasmic_df.tsv.gz"
    File cluster_coverage_txt_gz = "cluster/final/cluster.coverage.txt.gz"
    File cluster_depthTable_txt = "cluster/final/cluster.depthTable.txt"
    File cluster_rds = "cluster/final/cluster.rds"
    File cluster_signac_rds = "cluster/final/cluster.signac.rds"
    File cluster_variant_stats_tsv_gz = "cluster/final/cluster.variant_stats.tsv.gz"
    File cluster_vmr_strand_plot_png = "cluster/final/cluster.vmr_strand_plot.png"
    File passingBarcodes_tsv = "cluster/final/passingBarcodes.tsv"

  }
}


task plot_scmtah {

}