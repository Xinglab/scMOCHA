
workflow SCMTAH {
  # version of this pipeline
  String version = "CellRanger v7.0.1"

  String output_id
  String fastqs
  String sample_id
  String transcriptome
  String output_dir
  File rCRS

  Int? expect_cells

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
        memory = memory,
        boot_disk_size_gb = boot_disk_size_gb,
        disk_space = disk_space,
        cpu = cpu
  }

  call cell_cluster_annotation {
    input:
      h5file = cellranger_count.filtered_feature_bc_matrix,
      memory = memory,
      boot_disk_size_gb = boot_disk_size_gb,
      disk_space = disk_space,
      cpu = cpu
  }

  # call call_variant_on_single_cell_level {
  #   input:
  #   possorted_genome_bam = cellranger_count.sorted_bam,
  #   gzipped_barcodes = cellranger_count.barcodes,
  #   rCRS = rCRS,
  #   memory = memory,
  #   boot_disk_size_gb = boot_disk_size_gb,
  #   disk_space = disk_space,
  #   cpu = cpu
  # }

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
      File features = cellranger_count.features
      File matrix = cellranger_count.matrix


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
    Int? expect_cells

    String memory
    Int boot_disk_size_gb
    String disk_space
    Int cpu

    command {

      cellranger count \
        --id=${output_id} \
        --fastqs=${fastqs} \
        --sample=${sample_id} \
        --transcriptome=${transcriptome} \
        --nosecondary \
        --disable-ui \
        --localcores ${cpu}
    }

    output {
        File filtered_feature_bc_matrix = "${output_id}/outs/filtered_feature_bc_matrix.h5"
        File metrics_summary = "${output_id}/outs/metrics_summary.csv"
        File sorted_bam = "${output_id}/outs/possorted_genome_bam.bam"
        File sorted_bam_index = "${output_id}/outs/possorted_genome_bam.bam.bai"
        File web_summary = "${output_id}/outs/web_summary.html"
        File barcodes = "${output_id}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
        File features = "${output_id}/outs/filtered_feature_bc_matrix/features.tsv.gz"
        File matrix = "${output_id}/outs/filtered_feature_bc_matrix/matrix.mtx.gz"
    }
}

task cell_cluster_annotation {
  File h5file
  String refname = "pbmcref"
  String celllevel = "celltype.l1"

  String memory
  Int boot_disk_size_gb
  String disk_space
  Int cpu

  command {
    module load R/4.1.0
    Rscript /home/liuc9/github/scRNAseq-MitoVariant/bin/azimuth.R ${h5file} ${refname} ${celllevel}
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
  }
}

task call_variant_on_single_cell_level {
  File possorted_genome_bam
  File gzipped_barcodes
  File rCRS
  # Set 'barcodes' to the uncompressed barcode file
  File barcodes = sub(".gz$", "", gzipped_barcodes)


  String memory
  Int boot_disk_size_gb
  String disk_space
  Int cpu

  command {
    gunzip -c ${gzipped_barcodes} > ${barcodes}

    mgatk tenx -i ${possorted_genome_bam} \
      -n sc \
      -c ${cpu} -ub UB  -bt CB \
      -b ${barcodes} \
      --mito-genome ${rCRS}

    tar czf mgatk_single_cell_level.tar.gz "final"
  }
  output {
    File mgatk_single_cell_level = "mgatk_single_cell_level.tar.gz"
    File cell_heteroplasmic_df =  "final/sc.cell_heteroplasmic_df.tsv.gz"
    File coverage =  "final/sc.coverage.txt.gz"
  }
}

task call_variant_on_cell_cluster_level {
  File possorted_genome_bam
  File barcode_cluster
  File barcode_bulk
  File rCRS

  String memory
  Int boot_disk_size_gb
  String disk_space
  Int cpu

  command {

    samtools view -hb ${possorted_genome_bam} MT > MT.bam
    samtools index MT.bam

    # cell cluster

    sinto addtags \
      -b MT.bam \
      -f ${barcode_cluster} \
      -o MT_cluster.bam \
      -p ${cpu}

    samtools index MT_cluster.bam

    mgatk bcall -i MT_cluster.bam \
    -o mgatk_cluster \
    -n mgatk_cluster \
    -c ${cpu} -bt CJ \
    --mito-genome ${rCRS} \
    --keep-temp-files

    python /home/liuc9/github/scRNAseq-MitoVariant/bin/variant_calling.py mgatk_cluster/final/ mgatk_cluster 16569 10 MT


    # cell bulk

    sinto addtags \
      -b MT.bam \
      -f ${barcode_bulk} \
      -o MT_bulk.bam \
      -p ${cpu}

    samtools index MT_bulk.bam


    mgatk bcall -i MT_bulk.bam \
      -o mgatk_bulk \
      -n mgatk_bulk \
      -c ${cpu} -bt CJ \
      --mito-genome ${rCRS} \
      --keep-temp-files

    python /home/liuc9/github/scRNAseq-MitoVariant/bin/variant_calling.py mgatk_bulk/final/ mgatk_bulk 16569 10 MT
  }

}

task heteroplasmy_heatmap {
  File hetero_file
  File coverage_file
  File cluster_umap_file

  command {
    Rscript /home/liuc9/github/scRNAseq-MitoVariant/bin/heteroplasmy_heatmap.R ${hetero_file} ${coverage_file} ${cluster_umap_file}
  }
}