
workflow scMOCHA {
  # version of this pipeline
  String version = "CellRanger v7.0.1"

  # MitoScape
  Boolean use_mitoscape = false

  # Cell ranger inpunts
  String output_id
  String fastqs
  String sample_id
  String transcriptome = "/home/liuc9/data/refdata/mgatk_index/Human"
  File rCRS = "/home/liuc9/github/scMOCHA/fasta/rCRS.MT.fasta"

  File mt_exons_df = "/home/liuc9/github/scMOCHA/fasta/mt_exons.df.rds.gz"
  File mt_features_gmoviz = "/home/liuc9/github/scMOCHA/fasta/mt_features.grange.gmoviz.rds.gz"


  String output_dir

  # mgatk inputs
  String chrM = "MT"
  Int low_coverage_threshold = 10

  # cell_cluster_annotation inputs
  Int npcs = 10
  Float reso = 0.1
  String cellrefname
  String celllevel
  Int nFeature_RNA_min = 200
  Int nFeature_RNA_max = 8000
  Float percent_mt_max = 75
  Float percent_ribo_max = 50
  Float percent_Lagest_Gene_max = 50


  # Runtime attributes
  String memory = "50 GB"
  Int boot_disk_size_gb = 12
  String disk_space = "50"
  Int cpu = 10
  Boolean use_ssd = false

  # docker image
  String scmocha_version = "latest"
  String docker = "chunjiesamliu/scmocha"
  String partition = "defq"
  String account = "liuc9"
  File IMAGE = "/scr1/users/liuc9/sif/scmocha_latest.sif"

  File perlscript = "/home/liuc9/github/scMOCHA/bin/get_variants_info.pl"
  File jar_path = "/scr1/users/liuc9/tools/haplogrep3"  # /opt/haplogrep3/haplogrep3.jar
  File sqlite_path = "/mnt/isilon/xing_lab/liuc9/refdata/mitomaster/mitomap_sqlite_20230525.sqlite3"

  String bindir = "/home/liuc9/github/scMOCHA/bin"
  String conda_root = "/home/liuc9/tools/anaconda3"
  String conda_env = "scmocha"


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
        mt_exons_df = mt_exons_df,
        memory = memory,
        boot_disk_size_gb = boot_disk_size_gb,
        disk_space = disk_space,
        cpu = cpu,
        scmocha_version = scmocha_version,
        docker = docker,
        partition = partition,
        account = account,
        IMAGE = IMAGE,
        bindir = bindir,
        conda_root = conda_root,
        conda_env = conda_env

  }

  call cell_cluster_annotation {
    input:
      h5file = cellranger_count.filtered_feature_bc_matrix,
      mt_bam = cellranger_count.mt_bam,
      mt_bam_index = cellranger_count.mt_bam_index,
      npcs = npcs,
      reso = reso,
      refname = cellrefname,
      celllevel = celllevel,
      nFeature_RNA_min = nFeature_RNA_min,
      nFeature_RNA_max = nFeature_RNA_max,
      percent_mt_max = percent_mt_max,
      percent_ribo_max = percent_ribo_max,
      percent_Lagest_Gene_max = percent_Lagest_Gene_max,
      mt_rcrs_fasta = rCRS,
      mt_exons_df = mt_exons_df,
      mt_features_gmoviz = mt_features_gmoviz,
      memory = memory,
      boot_disk_size_gb = boot_disk_size_gb,
      disk_space = disk_space,
      cpu = cpu,
      scmocha_version = scmocha_version,
      docker = docker,
      partition = partition,
      account = account,
      IMAGE = IMAGE,
      bindir = bindir,
      conda_root = conda_root,
      conda_env = conda_env
  }

  call call_mt_variants {
    input:
      sorted_bam = cellranger_count.sorted_bam,
      sorted_bam_index = cellranger_count.sorted_bam_index,
      barcodes_tsv = cell_cluster_annotation.barcode_cell,
      mt_bam = cellranger_count.mt_bam,
      mt_bam_index = cellranger_count.mt_bam_index,
      mt_cluster_bam = cell_cluster_annotation.mt_cluster_bam,
      mt_cluster_bam_index = cell_cluster_annotation.mt_cluster_bam_index,
      mt_bulk_bam = cell_cluster_annotation.mt_bulk_bam,
      mt_bulk_bam_index = cell_cluster_annotation.mt_bulk_bam_index,
      chrM = chrM,
      rCRS = rCRS,
      low_coverage_threshold = low_coverage_threshold,
      memory = memory,
      boot_disk_size_gb = boot_disk_size_gb,
      disk_space = disk_space,
      cpu = cpu,
      scmocha_version = scmocha_version,
      docker = docker,
      partition = partition,
      account = account,
      IMAGE = IMAGE,
      bindir = bindir,
      conda_root = conda_root,
      conda_env = conda_env
  }

  call plot_scMOCHA {
    input:
      barcode_cluster_file = cell_cluster_annotation.barcode_cluster,
      cell_hetero_file = call_mt_variants.cell_cell_heteroplasmic_df_tsv_gz,
      cell_coverage_file = call_mt_variants.cell_coverage_txt_gz,
      cluster_hetero_file = call_mt_variants.cluster_cell_heteroplasmic_df_tsv_gz,
      cluster_coverage_file = call_mt_variants.cluster_coverage_txt_gz,
      cell_hetero_raw_file = call_mt_variants.cell_cell_heteroplasmic_df_raw_tsv_gz,
      perlscript = perlscript,
      jar_path = jar_path,
      sqlite_path = sqlite_path,
      memory = memory,
      boot_disk_size_gb = boot_disk_size_gb,
      disk_space = disk_space,
      cpu = cpu,
      scmocha_version = scmocha_version,
      docker = docker,
      partition = partition,
      account = account,
      IMAGE = IMAGE,
      bindir = bindir,
      conda_root = conda_root,
      conda_env = conda_env
  }

  call gather_outputfiles {
    input:
      output_dir = output_dir,
      # plot_scmth
      scMOCHA_rda = plot_scMOCHA.scMOCHA_rda,
      cell_af_heatmap = plot_scMOCHA.cell_af_heatmap,
      cell_depth_heatmap = plot_scMOCHA.cell_depth_heatmap,
      cluster_af_heatmap = plot_scMOCHA.cluster_af_heatmap,
      cluster_depth_heatmap = plot_scMOCHA.cluster_depth_heatmap,
      cluster_cell_af_heatmap = plot_scMOCHA.cluster_cell_af_heatmap,
      cluster_cell_depth_heatmap = plot_scMOCHA.cluster_cell_depth_heatmap,
      cell_variant_annotation_tsv = plot_scMOCHA.cell_variant_annotation_tsv,
      cell_variant_annotation_xlsx = plot_scMOCHA.cell_variant_annotation_xlsx,
      # call_mt_variants
      # cell level
      cell_A_txt_gz = call_mt_variants.cell_A_txt_gz,
      cell_C_txt_gz = call_mt_variants.cell_C_txt_gz,
      cell_G_txt_gz = call_mt_variants.cell_G_txt_gz,
      cell_T_txt_gz = call_mt_variants.cell_T_txt_gz,
      cell_cell_heteroplasmic_df_tsv_gz = call_mt_variants.cell_cell_heteroplasmic_df_tsv_gz,
      cell_cell_heteroplasmic_df_raw_tsv_gz = call_mt_variants.cell_cell_heteroplasmic_df_raw_tsv_gz,
      cell_coverage_txt_gz = call_mt_variants.cell_coverage_txt_gz,
      cell_depthTable_txt = call_mt_variants.cell_depthTable_txt,
      cell_rds = call_mt_variants.cell_rds,
      cell_signac_rds = call_mt_variants.cell_signac_rds,
      cell_variant_stats_tsv_gz = call_mt_variants.cell_variant_stats_tsv_gz,
      cell_vmr_strand_plot_png = call_mt_variants.cell_vmr_strand_plot_png,
      # cluster level
      barcodeQuants_tsv = call_mt_variants.barcodeQuants_tsv,
      cluster_A_txt_gz = call_mt_variants.cluster_A_txt_gz,
      cluster_C_txt_gz = call_mt_variants.cluster_C_txt_gz,
      cluster_G_txt_gz = call_mt_variants.cluster_G_txt_gz,
      cluster_T_txt_gz = call_mt_variants.cluster_T_txt_gz,
      cluster_cell_heteroplasmic_df_tsv_gz = call_mt_variants.cluster_cell_heteroplasmic_df_tsv_gz,
      cluster_coverage_txt_gz = call_mt_variants.cluster_coverage_txt_gz,
      cluster_depthTable_txt = call_mt_variants.cluster_depthTable_txt,
      cluster_rds = call_mt_variants.cluster_rds,
      cluster_signac_rds = call_mt_variants.cluster_signac_rds,
      cluster_variant_stats_tsv_gz = call_mt_variants.cluster_variant_stats_tsv_gz,
      cluster_vmr_strand_plot_png = call_mt_variants.cluster_vmr_strand_plot_png,
      passingBarcodes_tsv = call_mt_variants.passingBarcodes_tsv,
      # cell_cluster_annotation
      azimuth_rda = cell_cluster_annotation.azimuth_rda,
      barcode_cluster = cell_cluster_annotation.barcode_cluster,
      barcode_bulk = cell_cluster_annotation.barcode_bulk,
      celltype_ratio = cell_cluster_annotation.celltype_ratio,
      plot_metrics = cell_cluster_annotation.plot_metrics,
      plot_pie_celltype = cell_cluster_annotation.plot_pie_celltype,
      plot_qc = cell_cluster_annotation.plot_qc,
      plot_umap = cell_cluster_annotation.plot_umap,
      qc_cell_stats = cell_cluster_annotation.qc_cell_stats,
      sc_azimuth_rds_gz = cell_cluster_annotation.sc_azimuth_rds_gz,
      mt_cluster_bam = cell_cluster_annotation.mt_cluster_bam,
      mt_cluster_bam_index = cell_cluster_annotation.mt_cluster_bam_index,
      plot_mt_cluster_depth = cell_cluster_annotation.plot_mt_cluster_depth,
      # cellranger_count
      filtered_feature_bc_matrix = cellranger_count.filtered_feature_bc_matrix,
      metrics_summary = cellranger_count.metrics_summary,
      web_summary = cellranger_count.web_summary,
      barcodes = cellranger_count.barcodes,
      barcodes_tsv = cellranger_count.barcodes_tsv,
      features = cellranger_count.features,
      matrix = cellranger_count.matrix,
      mt_depth = cellranger_count.mt_depth,
      mt_depth_plot = cellranger_count.mt_depth_plot,
      mt_bam = cellranger_count.mt_bam,
      mt_bam_index = cellranger_count.mt_bam_index,
      memory = memory,
      boot_disk_size_gb = boot_disk_size_gb,
      disk_space = disk_space,
      cpu = cpu,
      scmocha_version = scmocha_version,
      docker = docker,
      partition = partition,
      account = account,
      IMAGE = IMAGE,
      bindir = bindir,
      conda_root = conda_root,
      conda_env = conda_env


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
      File mt_cluster_bam = cell_cluster_annotation.mt_cluster_bam
      File mt_cluster_bam_index = cell_cluster_annotation.mt_cluster_bam_index
      File mt_bulk_bam = cell_cluster_annotation.mt_bulk_bam
      File mt_bulk_bam_index = cell_cluster_annotation.mt_bulk_bam_index

      # call_mt_variants
      # cell level
      File cell_A_txt_gz = call_mt_variants.cell_A_txt_gz
      File cell_C_txt_gz = call_mt_variants.cell_C_txt_gz
      File cell_G_txt_gz = call_mt_variants.cell_G_txt_gz
      File cell_T_txt_gz = call_mt_variants.cell_T_txt_gz
      File cell_cell_heteroplasmic_df_tsv_gz = call_mt_variants.cell_cell_heteroplasmic_df_tsv_gz
      File cell_cell_heteroplasmic_df_raw_tsv_gz = call_mt_variants.cell_cell_heteroplasmic_df_raw_tsv_gz
      File cell_coverage_txt_gz = call_mt_variants.cell_coverage_txt_gz
      File cell_depthTable_txt = call_mt_variants.cell_depthTable_txt
      File cell_rds = call_mt_variants.cell_rds
      File cell_signac_rds = call_mt_variants.cell_signac_rds
      File cell_variant_stats_tsv_gz = call_mt_variants.cell_variant_stats_tsv_gz
      File cell_vmr_strand_plot_png = call_mt_variants.cell_vmr_strand_plot_png

      # cluster level
      File barcodeQuants_tsv = call_mt_variants.barcodeQuants_tsv
      File cluster_A_txt_gz = call_mt_variants.cluster_A_txt_gz
      File cluster_C_txt_gz = call_mt_variants.cluster_C_txt_gz
      File cluster_G_txt_gz = call_mt_variants.cluster_G_txt_gz
      File cluster_T_txt_gz = call_mt_variants.cluster_T_txt_gz
      File cluster_cell_heteroplasmic_df_tsv_gz = call_mt_variants.cluster_cell_heteroplasmic_df_tsv_gz
      File cluster_coverage_txt_gz = call_mt_variants.cluster_coverage_txt_gz
      File cluster_depthTable_txt = call_mt_variants.cluster_depthTable_txt
      File cluster_rds = call_mt_variants.cluster_rds
      File cluster_signac_rds = call_mt_variants.cluster_signac_rds
      File cluster_variant_stats_tsv_gz = call_mt_variants.cluster_variant_stats_tsv_gz
      File cluster_vmr_strand_plot_png = call_mt_variants.cluster_vmr_strand_plot_png
      File passingBarcodes_tsv = call_mt_variants.passingBarcodes_tsv

      # plot scMOCHA
      File scMOCHA_rda = plot_scMOCHA.scMOCHA_rda
      File cell_af_heatmap = plot_scMOCHA.cell_af_heatmap
      File cell_depth_heatmap = plot_scMOCHA.cell_depth_heatmap
      File cluster_af_heatmap = plot_scMOCHA.cluster_af_heatmap
      File cluster_depth_heatmap = plot_scMOCHA.cluster_depth_heatmap
      File cluster_cell_af_heatmap = plot_scMOCHA.cluster_cell_af_heatmap
      File cluster_cell_depth_heatmap = plot_scMOCHA.cluster_cell_depth_heatmap
      File cell_variant_annotation_tsv = plot_scMOCHA.cell_variant_annotation_tsv
      File cell_variant_annotation_xlsx = plot_scMOCHA.cell_variant_annotation_xlsx

      # gather_outputfiles
      File output_dir_tar_gz = gather_outputfiles.output_dir_tar_gz

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

    File mt_exons_df

    String memory
    Int boot_disk_size_gb
    String disk_space
    Int cpu
    String scmocha_version
    String docker
    String partition
    String account
    File IMAGE
    String bindir
    String conda_root
    String conda_env


    command {

      # module load R/4.1.0
      # module load R/4.2.3
      source ${conda_root}/etc/profile.d/conda.sh
      conda activate ${conda_env}

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

      # MT depth
      samtools depth -a -r ${chrM} --threads=${cpu} ${output_id}/outs/possorted_genome_bam.MT.bam > ${output_id}/outs/possorted_genome_bam.MT.depth

      # Depth plot
      ${bindir}/depth.R \
        -d ${output_id}/outs/possorted_genome_bam.MT.depth \
        -o ${output_id}/outs/possorted_genome_bam.MT.depth.pdf \
        -m ${mt_exons_df}

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

  Int npcs = 10
  Float reso = 0.1
  String refname
  String celllevel
  Int nFeature_RNA_min = 200
  Int nFeature_RNA_max = 8000
  Float percent_mt_max = 75
  Float percent_ribo_max = 50
  Float percent_Lagest_Gene_max = 50

  File mt_rcrs_fasta
  File mt_exons_df
  File mt_features_gmoviz

  String memory
  Int boot_disk_size_gb
  String disk_space
  Int cpu
  String scmocha_version
  String docker
  String partition
  String account
  File IMAGE
  String bindir
  String conda_root
  String conda_env


  command {
    # module load R/4.1.0
    # module load R/4.2.3
    source ${conda_root}/etc/profile.d/conda.sh
    conda activate ${conda_env}

    # cell cluster annotation
    ${bindir}/azimuth.R \
      -h5file ${h5file} \
      -npcs ${npcs} \
      -reso ${reso} \
      -refname_celllevel refname=${refname} celllevel=${celllevel} \
      -nFeature_RNA_min ${nFeature_RNA_min} \
      -nFeature_RNA_max ${nFeature_RNA_max} \
      -percent_mt_max ${percent_mt_max} \
      -percent_ribo_max ${percent_ribo_max} \
      -percent_Lagest_Gene_max ${percent_Lagest_Gene_max}


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

    # split MT bam by cluster
    bamtools split -in MT_cluster.bam -tag CJ

    # gmoviz plot of cluster coverage
    ${bindir}/depth_cluster_gmoviz.R ${mt_features_gmoviz} ${mt_rcrs_fasta} ${mt_exons_df}

  }
  output {
    File azimuth_rda = "azimuth.rda"
    File barcode_cell = "barcode_cell.tsv"
    File barcode_cluster = "barcode_cluster.tsv"
    File barcode_bulk = "barcode_bulk.tsv"
    File celltype_ratio = "celltype_ratio.tsv"
    File plot_metrics = "plot-metrics.pdf"
    File plot_pie_celltype = "plot-pie-celltype.pdf"
    File plot_qc = "plot-qc.pdf"
    File plot_umap = "plot-umap.pdf"
    File qc_cell_stats = "qc_cell_stats.xlsx"
    File sc_azimuth_rds_gz = "sc_azimuth.rds.gz"
    File mt_cluster_bam = "MT_cluster.bam"
    File mt_cluster_bam_index = "MT_cluster.bam.bai"
    File mt_bulk_bam = "MT_bulk.bam"
    File mt_bulk_bam_index = "MT_bulk.bam.bai"
    File plot_mt_cluster_depth = "plot-mt-cluster-depth.pdf"
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
  Int low_coverage_threshold


  String memory
  Int boot_disk_size_gb
  String disk_space
  Int cpu
  String scmocha_version
  String docker
  String partition
  String account
  File IMAGE
  String bindir
  String conda_root
  String conda_env


  command {

    # module load R/4.1.0
    # module load R/4.2.3
    source ${conda_root}/etc/profile.d/conda.sh
    conda activate ${conda_env}
    # call variants on single cell level
    mgatk tenx -i ${sorted_bam} \
      -o cell \
      -n cell \
      -g ${rCRS} \
      -c ${cpu} \
      -bt CB \
      -b ${barcodes_tsv} \
      -ub UB \
      --low-coverage-threshold ${low_coverage_threshold} \
      --snake-stdout \
      --keep-temp-files

    # the cell/final/ last "/" is important
    ${bindir}/variant_calling_cell_raw.py \
      cell/final/ \
      cell \
      16569 \
      ${low_coverage_threshold} \
      ${chrM}

    # call variants on cluster level
    mgatk bcall -i ${mt_cluster_bam} \
      -o cluster \
      -n cluster \
      -g ${rCRS} \
      -c ${cpu} \
      -bt CJ \
      --low-coverage-threshold ${low_coverage_threshold} \
      --snake-stdout \
      --keep-temp-files

    ${bindir}/variant_calling_cluster.py \
      cluster/final/ \
      cluster \
      16569 \
      ${low_coverage_threshold} \
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

task plot_scMOCHA {

  File barcode_cluster_file
  File cell_hetero_file
  File cell_coverage_file

  File cluster_hetero_file
  File cluster_coverage_file

  File cell_hetero_raw_file

  File perlscript
  File jar_path
  File sqlite_path

  String memory
  Int boot_disk_size_gb
  String disk_space
  Int cpu
  String scmocha_version
  String docker
  String partition
  String account
  File IMAGE
  String bindir
  String conda_root
  String conda_env


  command {
    # module load R/4.1.0
    # module load R/4.2.3
    source ${conda_root}/etc/profile.d/conda.sh
    conda activate ${conda_env}

    ${bindir}/scMOCHA.R \
      -b ${barcode_cluster_file} \
      -ceh ${cell_hetero_file} \
      -cec ${cell_coverage_file} \
      -clh ${cluster_hetero_file} \
      -clc ${cluster_coverage_file} \
      -chr ${cell_hetero_raw_file} \
      -p ${perlscript} \
      -j ${jar_path} \
      -s ${sqlite_path} \
      -conda_root ${conda_root} \
      -conda_env ${conda_env}
  }

  output {
    File scMOCHA_rda = "scMOCHA.rda"
    File cell_af_heatmap = "cell_af_heatmap.pdf"
    File cell_depth_heatmap = "cell_depth_heatmap.pdf"
    File cluster_af_heatmap = "cluster_af_heatmap.pdf"
    File cluster_depth_heatmap = "cluster_depth_heatmap.pdf"
    File cluster_cell_af_heatmap = "cluster_cell_af_heatmap.pdf"
    File cluster_cell_depth_heatmap = "cluster_cell_depth_heatmap.pdf"
    File cell_variant_annotation_tsv = "cell_variant_annotation.tsv"
    File cell_variant_annotation_xlsx = "cell_variant_annotation.xlsx"
  }

}

task gather_outputfiles {
  String output_dir

  # plot_scMOCHA
  File scMOCHA_rda
  File cell_af_heatmap
  File cell_depth_heatmap
  File cluster_af_heatmap
  File cluster_depth_heatmap
  File cluster_cell_af_heatmap
  File cluster_cell_depth_heatmap

  # call_mt_variants
  # cell level
  File cell_A_txt_gz
  File cell_C_txt_gz
  File cell_G_txt_gz
  File cell_T_txt_gz
  File cell_cell_heteroplasmic_df_tsv_gz
  File cell_cell_heteroplasmic_df_raw_tsv_gz
  File cell_coverage_txt_gz
  File cell_depthTable_txt
  File cell_rds
  File cell_signac_rds
  File cell_variant_stats_tsv_gz
  File cell_vmr_strand_plot_png
  # cluster level
  File barcodeQuants_tsv
  File cluster_A_txt_gz
  File cluster_C_txt_gz
  File cluster_G_txt_gz
  File cluster_T_txt_gz
  File cluster_cell_heteroplasmic_df_tsv_gz
  File cluster_coverage_txt_gz
  File cluster_depthTable_txt
  File cluster_rds
  File cluster_signac_rds
  File cluster_variant_stats_tsv_gz
  File cluster_vmr_strand_plot_png
  File passingBarcodes_tsv

  # cell_cluster_annotation
  File azimuth_rda
  File barcode_cluster
  File barcode_bulk
  File celltype_ratio
  File plot_metrics
  File plot_pie_celltype
  File plot_qc
  File plot_umap
  File qc_cell_stats
  File sc_azimuth_rds_gz
  File mt_cluster_bam
  File mt_cluster_bam_index
  File plot_mt_cluster_depth

  # cellranger_count
  File filtered_feature_bc_matrix
  File metrics_summary
  File web_summary
  File barcodes_tsv
  File barcodes
  File features
  File matrix
  File mt_depth
  File mt_depth_plot
  File mt_bam
  File mt_bam_index
  File cell_variant_annotation_tsv
  File cell_variant_annotation_xlsx

  String memory
  Int boot_disk_size_gb
  String disk_space
  Int cpu
  String scmocha_version
  String docker
  String partition
  String account
  File IMAGE
  String bindir
  String conda_root
  String conda_env

  command {
    mkdir -p ${output_dir}
    # plot_scMOCHA
    cp ${scMOCHA_rda} ${output_dir}
    cp ${cell_af_heatmap} ${output_dir}
    cp ${cell_depth_heatmap} ${output_dir}
    cp ${cluster_af_heatmap} ${output_dir}
    cp ${cluster_depth_heatmap} ${output_dir}
    cp ${cluster_cell_af_heatmap} ${output_dir}
    cp ${cluster_cell_depth_heatmap} ${output_dir}

    # call_mt_variants
    # cell level
    cp ${cell_A_txt_gz} ${output_dir}
    cp ${cell_C_txt_gz} ${output_dir}
    cp ${cell_G_txt_gz} ${output_dir}
    cp ${cell_T_txt_gz} ${output_dir}
    cp ${cell_cell_heteroplasmic_df_tsv_gz} ${output_dir}
    cp ${cell_cell_heteroplasmic_df_raw_tsv_gz} ${output_dir}
    cp ${cell_coverage_txt_gz} ${output_dir}
    cp ${cell_depthTable_txt} ${output_dir}
    cp ${cell_rds} ${output_dir}
    cp ${cell_signac_rds} ${output_dir}
    cp ${cell_variant_stats_tsv_gz} ${output_dir}
    cp ${cell_vmr_strand_plot_png} ${output_dir}
    # cluster level
    cp ${barcodeQuants_tsv} ${output_dir}
    cp ${cluster_A_txt_gz} ${output_dir}
    cp ${cluster_C_txt_gz} ${output_dir}
    cp ${cluster_G_txt_gz} ${output_dir}
    cp ${cluster_T_txt_gz} ${output_dir}
    cp ${cluster_cell_heteroplasmic_df_tsv_gz} ${output_dir}
    cp ${cluster_coverage_txt_gz} ${output_dir}
    cp ${cluster_depthTable_txt} ${output_dir}
    cp ${cluster_rds} ${output_dir}
    cp ${cluster_signac_rds} ${output_dir}
    cp ${cluster_variant_stats_tsv_gz} ${output_dir}
    cp ${cluster_vmr_strand_plot_png} ${output_dir}
    cp ${passingBarcodes_tsv} ${output_dir}

    # cell_cluster_annotation
    cp ${azimuth_rda} ${output_dir}
    cp ${barcode_cluster} ${output_dir}
    cp ${barcode_bulk} ${output_dir}
    cp ${celltype_ratio} ${output_dir}
    cp ${plot_metrics} ${output_dir}
    cp ${plot_pie_celltype} ${output_dir}
    cp ${plot_qc} ${output_dir}
    cp ${plot_umap} ${output_dir}
    cp ${qc_cell_stats} ${output_dir}
    cp ${sc_azimuth_rds_gz} ${output_dir}
    cp ${mt_cluster_bam} ${output_dir}
    cp ${mt_cluster_bam_index} ${output_dir}
    cp ${plot_mt_cluster_depth} ${output_dir}

    # cellranger_count
    cp ${filtered_feature_bc_matrix} ${output_dir}
    cp ${metrics_summary} ${output_dir}
    cp ${web_summary} ${output_dir}
    cp ${barcodes_tsv} ${output_dir}
    cp ${barcodes} ${output_dir}
    cp ${features} ${output_dir}
    cp ${matrix} ${output_dir}
    cp ${mt_depth} ${output_dir}
    cp ${mt_depth_plot} ${output_dir}
    cp ${mt_bam} ${output_dir}
    cp ${mt_bam_index} ${output_dir}
    cp ${cell_variant_annotation_tsv} ${output_dir}
    cp ${cell_variant_annotation_xlsx} ${output_dir}


    tar -czf ${output_dir}.tar.gz ${output_dir}
  }

  output {
    File output_dir_tar_gz = "${output_dir}.tar.gz"
  }


}