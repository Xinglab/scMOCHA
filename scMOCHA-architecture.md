# scMOCHA Architecture — Figure Description for Nature Review Genetics-Style Illustration

## What is scMOCHA?

**scMOCHA** (Single-Cell Mitochondrial Omics for Cellular Heteroplasmy Analysis) is an end-to-end computational pipeline for calling, annotating, and visualizing **mitochondrial DNA (mtDNA) variants** from **10x Genomics single-cell RNA-seq** data. It operates at both single-cell and cell-cluster resolution, enabling detection of **heteroplasmy** — the co-existence of wild-type and mutant mtDNA within individual cells.

The pipeline is implemented in **WDL (Workflow Description Language)** and orchestrated by Cromwell, enabling execution on local machines, HPC clusters (SLURM), and cloud platforms. It supports 10x Genomics chemistry versions v1 through v5.

### Key Innovation: Simultaneous Multi-Modal Calling
A core advantage of scMOCHA is its ability to **simultaneously extract three dimensions of single-cell information from a single standard scRNA-seq library**:
1. **Mitochondrial Variants & Heteroplasmy**: Calling somatic mitochondrial DNA mutations at both individual cell and cell-type cluster levels.
2. **Nuclear Gene Expression**: Capturing host nuclear gene expression profiles to determine overall transcriptomic state.
3. **Cell Type Identity**: Identifying cell types and subpopulations (via reference-based Azimuth, marker-based scType, or unsupervised Seurat clustering) using the nuclear gene expression data.

By calling variants, expression, and cell type simultaneously from the same sequence library, scMOCHA links genetic heteroplasmy directly to phenotypic cell states without requiring separate, complex multi-omic assays.

---

## Suggested Figure Layout

> **Overall layout**: A landscape-oriented figure divided into **4 horizontal panels (A–D)** arranged left-to-right, representing the pipeline flow. Each panel contains a shaded box grouping related steps. Arrows connect panels to show data flow. Below the main pipeline, a narrow **Panel E** spans the full width showing representative output visualizations.

```
┌─────────────────────────────────────────────────────────────────────────────────────┐
│  Panel A              Panel B                Panel C               Panel D          │
│  ┌──────────┐    ┌──────────────────┐    ┌───────────────┐    ┌─────────────────┐   │
│  │ INPUT &  │───▶│  CELL ANALYSIS   │───▶│  MT VARIANT   │───▶│  ANNOTATION &   │   │
│  │ALIGNMENT │    │  & MT READ       │    │  CALLING      │    │  VISUALIZATION  │   │
│  │          │    │  EXTRACTION      │    │               │    │                 │   │
│  └──────────┘    └──────────────────┘    └───────────────┘    └─────────────────┘   │
│                                                                                     │
│  ───────────────────────── Panel E: Representative Outputs ─────────────────────     │
│  ┌─────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐      │
│  │ QC Plot │ │  UMAP    │ │ Circular │ │ Heatmap  │ │ Violin   │ │ Feature  │      │
│  │         │ │          │ │ mtGenome │ │          │ │ Plot     │ │ Plot     │      │
│  └─────────┘ └──────────┘ └──────────┘ └──────────┘ └──────────┘ └──────────┘      │
└─────────────────────────────────────────────────────────────────────────────────────┘
```

---

## Panel A — Input & Alignment

### Title: **"Raw Data Processing"**
### Color theme: Cool blue / teal

This panel illustrates the starting point of the pipeline.

### Elements to draw:

1. **10x Genomics scRNA-seq FASTQ files** (icon: a stylized sequencing data file or a flow cell icon)
   - Show Read 1 (cell barcode + UMI, 26–28 bp), Read 2 (cDNA insert, 98 bp), and Index read (I7, 8 bp sample index)
   - Label: "10x scRNA-seq (v1/v2/v3/v3HT/v4/v5)"

2. **Chemistry Auto-Detection** (small icon: gear/settings)
   - `chemistry.R` reads FASTQ headers and determines the 10x chemistry version automatically
   - Outputs: chemistry version string (e.g., "SC3Pv3")

3. **Cell Ranger `count`** (icon: Cell Ranger logo or alignment arrows)
   - Input: FASTQ files + GRCh38 reference genome **with rCRS mitochondrial sequence**
   - Key detail: The standard mitochondrial reference is replaced with the **Revised Cambridge Reference Sequence (rCRS)** — the gold standard for human mtDNA coordinate system
   - Outputs: **BAM file** (aligned reads), **filtered barcode matrix** (gene × cell count matrix), **web summary** (QC HTML report)
   - Arrow pointing right: BAM goes to Panel B; count matrix goes to Panel B

### Key annotation text:
- "Reference genome: GRCh38 + rCRS (MT)"
- "Auto-detects 10x chemistry v1–v5"

---

## Panel B — Cell Analysis & Mitochondrial Read Extraction

### Title: **"Cell QC, Clustering & MT Read Extraction"**
### Color theme: Green / emerald

This panel has **two parallel tracks** flowing from Cell Ranger output:

### Track 1 (top): Cell Quality Control & Annotation

4. **Seurat QC Filtering** (icon: filter funnel)
   - Filters cells by:
     - `nFeature_RNA`: 200–8,000 genes per cell
     - `percent.mt` < 75% (mitochondrial content)
     - `percent.ribo` < 50% (ribosomal content)
     - `percent.largest_gene` < 50%
   - Small schematic: scatter plot of nCount_RNA vs nFeature_RNA with red threshold lines

5. **Cell Type Annotation** (icon: cell clusters)
   - **Azimuth** (primary method): Reference-based annotation for known tissue types (PBMC, bone marrow, brain cortex, lung, kidney, heart, adipose, etc.)
   - **scType** (alternative): Marker-based annotation using curated gene sets
   - **Standard Seurat clustering** (fallback): PCA → UMAP → FindNeighbors → FindClusters
   - Output: UMAP plot with color-coded cell types (e.g., CD4 T, CD8 T, Mono, B, NK, DC)
   - Output: **barcode-to-cluster mapping** (sent to both tracks)

### Track 2 (bottom): Mitochondrial Read Extraction

6. **samtools MT Read Extraction** (icon: DNA strand filtering)
   - `samtools view` extracts reads mapped to the MT chromosome from the BAM
   - Produces: **MT-only BAM file**

7. **sinto Cell-Barcode Splitting** (icon: branching arrows)
   - `sinto filterbarcodes` splits the MT BAM into per-cluster BAM files
   - Uses the barcode-to-cluster mapping from Track 1
   - Produces: Individual BAM files per cell-type cluster (e.g., `CD4_T.bam`, `Mono.bam`, `B.bam`)

### Visual connector:
- Dashed arrow from Track 1's barcode-cluster mapping to Track 2's sinto step, showing the integration point

### Key annotation text:
- "Dual-track: Cell annotation informs MT read partitioning"
- "Azimuth / scType / Seurat for automated cell typing"

---

## Panel C — Mitochondrial Variant Calling

### Title: **"Dual-Resolution Variant Calling"**
### Color theme: Warm orange / amber

This panel shows the **two parallel variant calling paths** — the defining feature of scMOCHA.

### Path 1 (top): Single-Cell Level

8. **mgatk (Single-Cell)** (icon: variant calling bars)
   - **mgatk** (mitochondrial genotyping toolkit) processes the cell-level MT BAM
   - Counts A, C, G, T nucleotides at every MT position (1–16,569) for every cell barcode
   - Tracks forward and reverse strand counts separately
   - Outputs: Sparse matrices (`.mtx` format) of per-base allele counts

9. **Cell-Level Variant Processing** (`variant_calling_cell_raw.py`)
   - Computes allele frequency: AF = alt_reads / total_reads per cell per position
   - Filters: minimum coverage ≥ 10 reads, strand support on both strands
   - Computes strand concordance: log₂(forward/reverse) to detect strand bias
   - Outputs: **Cell AF matrix** (variants × cells), **Cell depth matrix**, **variant list**

### Path 2 (bottom): Cluster Level

10. **mgatk (Cluster)** (icon: aggregated variant calling bars)
    - Runs mgatk separately on each cluster BAM file
    - Higher coverage per position (aggregated reads) → more sensitive variant detection
    - Outputs: Per-cluster allele count matrices

11. **Cluster-Level Variant Processing** (`variant_calling_cluster.py`)
    - Same algorithm as cell-level but with clusters as columns
    - Aggregated allele frequencies per cell-type cluster
    - Outputs: **Cluster AF matrix** (variants × clusters), **Cluster depth matrix**, **variant list**

### Visual element:
- Show the MT genome as a horizontal bar (16,569 bp) with ticks
- Overlay small bar charts representing allele frequency at variant positions
- One bar chart for single-cell (sparse, many cells) and one for cluster (dense, few clusters)

### Key annotation text:
- "Dual resolution: per-cell AND per-cluster heteroplasmy"
- "mgatk: position-level allele counting"
- "Strand concordance filter removes artifacts"

---

## Panel D — Annotation & Visualization

### Title: **"Biological Annotation & Integrative Visualization"**
### Color theme: Purple / magenta

This panel shows the final annotation and output generation steps.

### Annotation Layer

12. **MITOMAP Annotation** (icon: database with DNA helix)
    - Perl script (`get_variants_info.pl`) queries MITOMAP SQLite database
    - Retrieves for each variant:
      - **Locus**: Gene/region on MT genome (e.g., MT-ND1, MT-CO1, MT-RNR1, D-loop, tRNA genes)
      - **Nucleotide change type**: Transition (Ti) vs. Transversion (Tv)
      - **Conservation score**: 0–100 (cross-species evolutionary conservation)
      - **Disease association**: Known mitochondrial disease links (e.g., LHON, MELAS, hearing loss, diabetes)
      - **MITOMAP allele frequency**: Population frequency from MITOMAP database
      - **gnomAD MT frequency**: Population frequency from gnomAD mitochondrial database
    - Output: Annotated variant table (TSV + Excel)

13. **Haplogrep3 Classification** (icon: phylogenetic tree branch)
    - Java-based haplogroup classifier
    - Assigns mitochondrial haplogroup (e.g., H, T2b, J1c, U5a) to each cell/cluster based on detected variants
    - Uses the PhyloTree reference for haplogroup assignment
    - Output: Haplogroup labels per cell/cluster

### Visualization Layer

14. **Heteroplasmy Heatmap** (icon: heatmap grid) — `scMOCHA.R` + `ComplexHeatmap`
    - Rows: Variants (ordered by MT position)
    - Columns: Cells (ordered by cluster) or Clusters
    - Color: Allele frequency (white→red, 0→1)
    - Multi-layer annotation:
      - **Top**: Cell cluster color bar, MT%, log₁₀(total reads)
      - **Left**: MITOMAP freq, gnomAD freq, haplogroup
      - **Right**: Conservation, nucleotide change, locus, disease

15. **Violin + Beeswarm Plots** — Per-variant allele frequency distribution across cell types
    - Each facet = one variant
    - X-axis = cell types
    - Y-axis = allele frequency (0–1)
    - Points colored by AF (gradient: light → dark)

16. **Depth Coverage Plots** — `depth.R`
    - Area plots showing read depth along the MT genome (positions 1–16,569)
    - Faceted by cell type
    - Bottom track: Gene annotation (rRNAs in dark, tRNAs as diamonds, protein-coding genes in blue)

17. **Circular MT Genome Plot** — `depth_cluster_gmoviz.R`
    - Circular ideogram of the 16,569 bp mitochondrial genome
    - Outer ring: Gene/feature annotations (color-coded by biotype)
    - Inner rings: Coverage depth per cluster (color-coded tracks)

18. **UMAP Feature Plots** — `variant_feature_plot.R`
    - Standard UMAP embedding with cells colored by allele frequency of specific variants
    - Shows spatial patterns of heteroplasmy across cell populations

19. **Variant Count Plots** — `variant_count_plot.R`
    - Bar plots of heteroplasmic variant count per cell, grouped by cell type
    - Shows variant burden differences across cell populations

### Key annotation text:
- "Comprehensive annotation: MITOMAP + gnomAD + haplogroup"
- "6 output visualization types"

---

## Panel E — Representative Outputs (Bottom Strip)

A narrow horizontal strip showing **miniature versions** of the 6 main output types:

| Output | Description | Visual Style |
|--------|-------------|--------------|
| **QC Plot** | Scatter plot: nCount_RNA vs nFeature_RNA, colored by %MT; histogram of %MT distribution | Scatter + histogram |
| **UMAP** | Cell type clusters on 2D UMAP embedding, color-coded | Dot cloud with labeled clusters |
| **Circular MT Genome** | gmoviz circular plot with gene annotations and depth tracks | Circular ideogram |
| **Heteroplasmy Heatmap** | ComplexHeatmap with multi-layer annotations | Dense heatmap with sidebars |
| **Violin Plot** | Allele frequency distributions per variant per cell type | Faceted violin/beeswarm |
| **Feature Plot** | UMAP colored by variant allele frequency | Gradient-colored dot cloud |

---

## Complete Data Flow Summary

```
  ┌─────────────────────────────────────────────────────────────────────────────────────────────┐
  │                         scMOCHA Pipeline Architecture                                       │
  │                                                                                             │
  │  ┌───────────┐     ┌──────────────┐                                                        │
  │  │  FASTQ    │────▶│  Chemistry   │                                                        │
  │  │  Files    │     │  Detection   │                                                        │
  │  │ (10x     │     │  (chemistry.R)│                                                        │
  │  │  v1-v5)  │     └──────┬───────┘                                                        │
  │  └───────────┘           │                                                                 │
  │                          ▼                                                                 │
  │              ┌──────────────────────┐                                                      │
  │              │    Cell Ranger       │                                                      │
  │              │  (Alignment + Count) │                                                      │
  │              │  GRCh38 + rCRS ref   │                                                      │
  │              └────────┬─────────────┘                                                      │
  │                       │                                                                    │
  │          ┌────────────┼────────────────┐                                                   │
  │          │            │                │                                                    │
  │          ▼            │                ▼                                                    │
  │  ┌──────────────┐    │    ┌──────────────────────┐                                        │
  │  │  Cell QC &   │    │    │  MT Read Extraction   │                                        │
  │  │  Clustering  │    │    │  (samtools + sinto)    │                                        │
  │  │  (Seurat +   │    │    │                        │                                        │
  │  │   Azimuth)   │    │    │  ┌────────┐ ┌────────┐│                                        │
  │  │              │─────────▶  │Cell BAM│ │Cluster ││                                        │
  │  │ ┌──────────┐│    │    │  │        │ │BAMs    ││                                        │
  │  │ │ Barcode- ││    │    │  └───┬────┘ └───┬────┘│                                        │
  │  │ │ Cluster  ││    │    └──────┼──────────┼─────┘                                        │
  │  │ │ Mapping  ││    │           │          │                                               │
  │  │ └──────────┘│    │           ▼          ▼                                               │
  │  └──────────────┘    │   ┌──────────┐ ┌──────────┐                                        │
  │          │           │   │  mgatk   │ │  mgatk   │                                        │
  │          │           │   │ (Cell)   │ │(Cluster) │                                        │
  │          │           │   └────┬─────┘ └────┬─────┘                                        │
  │          │           │        │            │                                               │
  │          │           │        ▼            ▼                                               │
  │          │           │   ┌──────────┐ ┌──────────┐                                        │
  │          │           │   │  Cell    │ │ Cluster  │                                        │
  │          │           │   │ Variant  │ │ Variant  │                                        │
  │          │           │   │ Calling  │ │ Calling  │                                        │
  │          │           │   │ (Python) │ │ (Python) │                                        │
  │          │           │   └────┬─────┘ └────┬─────┘                                        │
  │          │           │        │            │                                               │
  │          │           │        └──────┬─────┘                                               │
  │          │           │               │                                                     │
  │          │           │        ┌──────┴──────┐                                              │
  │          │           │        ▼             ▼                                              │
  │          │           │  ┌──────────┐ ┌────────────┐                                       │
  │          │           │  │ MITOMAP  │ │ Haplogrep3 │                                       │
  │          │           │  │ Annot.   │ │ Haplogroup │                                       │
  │          │           │  │ (Perl)   │ │ (Java)     │                                       │
  │          │           │  └────┬─────┘ └─────┬──────┘                                       │
  │          │           │       │             │                                               │
  │          │           │       └──────┬──────┘                                               │
  │          │           │              │                                                      │
  │          ▼           │              ▼                                                      │
  │  ┌───────────────────────────────────────────┐                                            │
  │  │         Integrative Visualization          │                                            │
  │  │  (scMOCHA.R — R/ComplexHeatmap/ggplot2)    │                                            │
  │  │                                             │                                            │
  │  │  ┌────────┐ ┌────────┐ ┌────────┐          │                                            │
  │  │  │Heatmap │ │Violin  │ │ UMAP   │          │                                            │
  │  │  │(AF)    │ │(AF by  │ │Feature │          │                                            │
  │  │  │        │ │celltype│ │ Plot   │          │                                            │
  │  │  └────────┘ └────────┘ └────────┘          │                                            │
  │  │  ┌────────┐ ┌────────┐ ┌────────┐          │                                            │
  │  │  │Depth   │ │Circular│ │Variant │          │                                            │
  │  │  │Coverage│ │mtGenome│ │ Count  │          │                                            │
  │  │  │Plot    │ │(gmoviz)│ │Barplot │          │                                            │
  │  │  └────────┘ └────────┘ └────────┘          │                                            │
  │  └────────────────────────────────────────────┘                                            │
  └────────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## Style Guidance for Nature Review Genetics Figure

### General Style
- **Clean, vector-based illustration** — no pixelation, no screenshots
- **Muted, harmonious color palette** — use Nature Reviews' signature style:
  - Soft blues, teals, greens, ambers, purples
  - Avoid harsh primary colors
  - Use white/light grey backgrounds with subtle shadows for depth
- **Consistent typography**: Sans-serif (Helvetica or Arial), 8–10 pt for labels
- **Rounded-corner boxes** for computational steps
- **Thin, clean arrows** with arrowheads for data flow (1–2 pt weight)
- **Icons/pictograms** for biological elements (cells, DNA, mitochondria)
- **Panel labels**: Bold letters A, B, C, D, E in top-left corners

### Color Coding Suggestions
| Element | Color | Hex |
|---------|-------|-----|
| Panel A (Input/Alignment) | Steel blue | #4A90B8 |
| Panel B (Cell Analysis) | Sea green | #3CB371 |
| Panel C (Variant Calling) | Warm amber | #E8A735 |
| Panel D (Annotation/Vis) | Soft purple | #8B6BB5 |
| Reference databases | Pale sage | #B8D4B8 |
| Software tools | White boxes with colored borders | |
| Data files/matrices | Rounded grey pills | #D0D0D0 |
| Arrows (data flow) | Dark grey | #555555 |

### Iconography
- **Mitochondrion**: Oval with cristae inner membrane folds (use standard biology textbook style)
- **Cell clusters**: Overlapping colored circles (each color = a cell type)
- **Heatmap mini-icon**: Small grid with gradient fill
- **UMAP mini-icon**: Cloud of colored dots
- **DNA helix**: For variant calling section
- **Database cylinder**: For MITOMAP/gnomAD

### Key Design Elements from Nature Review Genetics Style
1. **Numbered steps** in circles (①②③...) connected by arrows
2. **Inset boxes** with tool names and brief descriptions
3. **Data transformation callouts**: Small diagrams showing matrix shapes (e.g., "cells × genes" → "variants × cells")
4. **Legend** in bottom-right with symbol definitions
5. **Scale-appropriate detail**: Main flow is readable at journal column width; details visible at full page width

---

## Detailed Text for Each Step (for Figure Labels)

### Step labels (short, for inside boxes):
1. **Chemistry Detection** — Auto-detect 10x version (v1–v5)
2. **Cell Ranger** — Alignment to GRCh38 + rCRS
3. **Cell QC** — Filter by nFeature, %MT, %ribo
4. **Cell Annotation** — Azimuth / scType / Seurat
5. **MT Read Extraction** — samtools + sinto
6. **mgatk (Cell)** — Per-cell allele counting
7. **mgatk (Cluster)** — Per-cluster allele counting
8. **Cell Variant Calling** — AF, depth, strand concordance
9. **Cluster Variant Calling** — Aggregated AF per cluster
10. **MITOMAP Annotation** — Locus, conservation, disease
11. **Haplogrep3** — Haplogroup classification
12. **Visualization** — Heatmaps, violin, UMAP, depth, circular

### Step descriptions (longer, for figure legend or callout boxes):
1. Reads FASTQ headers to determine barcode length and chemistry version
2. Aligns reads to human genome with rCRS mitochondrial reference; generates gene expression count matrix
3. Filters low-quality cells by gene count, mitochondrial content, ribosomal content, and dominant gene percentage
4. Automated cell type annotation using reference-based (Azimuth) or marker-based (scType) methods
5. Extracts MT-mapped reads; splits by cell barcode into per-cluster BAM files
6. Counts A/C/G/T alleles at all 16,569 MT positions per individual cell (forward + reverse strand)
7. Same as above but aggregated per cell-type cluster for higher sensitivity
8. Computes allele frequencies, filters by coverage (≥10x) and strand concordance
9. Computes cluster-level allele frequencies with same quality filters
10. Queries MITOMAP database for disease associations, conservation, locus, nucleotide change type, and gnomAD frequencies
11. Classifies mitochondrial haplogroup from variant profile using PhyloTree reference
12. Generates 6 visualization types: heteroplasmy heatmap, violin plot, UMAP feature plot, depth coverage, circular genome, variant count

---

## Software & Database Summary Table

| Category | Tool/Database | Version | Purpose |
|----------|---------------|---------|---------|
| **Alignment** | Cell Ranger | v7.0.1 | scRNA-seq alignment + counting |
| **QC/Clustering** | Seurat | v4/v5 | Cell QC, clustering, UMAP |
| **Annotation** | Azimuth | latest | Reference-based cell typing |
| **BAM Processing** | samtools | ≥1.15 | MT read extraction, depth calculation |
| **BAM Splitting** | sinto | latest | Split BAM by cell barcode |
| **Variant Calling** | mgatk | latest | MT allele counting per cell |
| **Variant Calling** | Python (scipy, pandas) | 3.x | AF computation, filtering |
| **Annotation** | MITOMAP (SQLite) | latest | Disease association, conservation |
| **Annotation** | gnomAD (MT) | latest | Population allele frequencies |
| **Haplogroup** | Haplogrep3 | latest | MT haplogroup classification |
| **Visualization** | ComplexHeatmap | latest | Multi-annotation heatmaps |
| **Visualization** | gmoviz | latest | Circular MT genome plots |
| **Visualization** | ggplot2 | latest | Violin, scatter, bar plots |
| **Orchestration** | Cromwell/WDL | latest | Workflow execution engine |

---

## Key Scientific Concepts to Highlight in the Figure

1. **Heteroplasmy**: The co-existence of multiple mtDNA genotypes within a single cell. scMOCHA quantifies this as **allele frequency (AF)** at each variant position per cell.

2. **Dual-resolution analysis**: Variants are called at both:
   - **Single-cell level**: Reveals cell-to-cell heteroplasmy variation (lower coverage, more noise)
   - **Cluster level**: Aggregates reads across cells of the same type (higher coverage, more confident calls)

3. **rCRS standard**: The Revised Cambridge Reference Sequence is the internationally agreed-upon reference for human mtDNA, ensuring consistent variant nomenclature (e.g., "73A>G", "16519T>C").

4. **Strand concordance**: A quality metric that checks if a variant is supported by reads from both the forward and reverse strand, filtering out sequencing artifacts.

5. **The mitochondrial genome**: A circular, 16,569 bp genome encoding:
   - 13 protein-coding genes (respiratory chain subunits)
   - 22 tRNA genes
   - 2 rRNA genes (12S and 16S)
   - 1 non-coding control region (D-loop)

6. **Simultaneous Multi-Modal Profiling**: Highlight the direct, simultaneous recovery of mitochondrial variants, nuclear transcriptomes, and cell-type metadata from the same single-cell library, linking genotype to phenotype.

---

## Prompt Template for Image Generation (Gemini / GPT-image-2)

> Create a publication-quality scientific figure in the style of Nature Review Genetics, showing the scMOCHA pipeline architecture. The figure should be landscape-oriented with 4 main panels (A–D) arranged left-to-right, plus a bottom output strip (Panel E).
>
> **Key Visual Focus**: Highlight that scMOCHA **simultaneously calls mitochondrial variants, nuclear gene expression, and cell type information** from the same single-cell RNA-seq run.
>
> **Panel A "Input & Alignment"** (steel blue theme): Shows 10x Genomics scRNA-seq FASTQ files entering Chemistry Detection (auto-detect v1-v5), then Cell Ranger alignment to GRCh38+rCRS reference genome. Output: BAM file + filtered barcode matrix.
>
> **Panel B "Cell Analysis & MT Extraction"** (sea green theme): Two parallel tracks. Top track: Seurat QC filtering (by gene count, %MT, %ribo) → Azimuth/scType cell type annotation (nuclear gene expression analysis) → UMAP with labeled cell clusters. Bottom track: samtools extracts MT reads → sinto splits into per-cluster BAMs. A prominent label/badge highlights the **simultaneous call of nuclear expression and cell type mapping** linked to MT read separation. Dashed arrow connects cell annotation to BAM splitting.
>
> **Panel C "Dual-Resolution Variant Calling"** (warm amber theme): Two parallel paths. Top: mgatk (cell-level) → Python variant processing → cell AF matrix. Bottom: mgatk (cluster-level) → Python variant processing → cluster AF matrix. Show the MT genome (16,569 bp linear bar) with variant positions marked. Accentuate the simultaneous output of cell-level variants matched to cell-type annotations.
>
> **Panel D "Annotation & Visualization"** (soft purple theme): MITOMAP database annotation (locus, conservation, disease) + Haplogrep3 haplogroup classification → feeds into 6 output visualization types showing integrated profiles of MT variant allele frequencies, cell-type clusters, and gene expression metrics.
>
> **Panel E (bottom strip)**: Miniature representative outputs — QC scatter plot, UMAP, circular MT genome, heteroplasmy heatmap, violin plot, UMAP feature plot.
>
> Style: Clean vector illustration, muted harmonious colors, rounded-corner boxes, thin arrows, sans-serif typography, numbered steps in circles, Nature Reviews aesthetic. Include a small mitochondrion icon, cell cluster icons, and database cylinder icons. White background with subtle shadows.
