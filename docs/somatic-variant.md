# Workflow for Somatic mtDNA Variant Analysis

A primary objective of scMOCHA is to identify bona fide somatic mtDNA variants arising post-zygotically within single cells, while strictly isolating them from germline polymorphisms, haplogroup-defining variants, and technical sequencing artifacts. Because somatic mtDNA variants typically exhibit low cellular prevalence and low allele frequencies within an individual, they cannot be reliably identified from aggregated pseudobulk data alone. To resolve this, scMOCHA integrates single-cell, per-cluster, and global bulk-level evidence into a unified data structure, applying layered filtration windows to confidently capture sub-clonal mutations.

---

## Cell-Level Evidence Matrix

For every identified sample–variant pair, each single cell is cross-referenced and classified into one of four mutually exclusive evidence categories based on localized read support:

| Category         | Label        | Criteria                                                                                                                                                      |
| ---------------- | ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Variant-Positive | **Colorful** | ≥2 alternative-allele reads on both forward and reverse strands; total position depth ≥10 reads; cellular allele frequency > detection threshold (e.g., 0.05) |
| Reference-Only   | **Black**    | Sufficient coverage at target position (≥10 reads) with zero supporting evidence for alternative alleles                                                      |
| Low Coverage     | **Grey**     | Uninformative or low-confidence coverage (1–9 reads), yielding an ambiguous status                                                                            |
| No Coverage      | **White**    | No overlapping read alignments at the specific coordinate                                                                                                     |

---

## Pre-Filtering of Germline and Technical Artifacts

Prior to somatic nomination, all tracked variants are evaluated against three independent exclusion manifests to remove inherited and technical variants:

1. **Haplogroup Signatures** — Lineage-defining variants are extracted from the per-sample [HaploGrep3](https://haplogrep.i-med.ac.at/) profile based on the PhyloTree-mt 17.2 baseline and immediately purged from the somatic candidate pool.

2. **Common Population Variants** — Variants matching known inherited alleles are filtered using the [gnomAD v3.1](https://gnomad.broadinstitute.org/) mitochondrial release. Loci with a population homoplasmic allele frequency (AF<sub>homo</sub>) exceeding a permissive cutoff, or those found as homoplasmies in ≥5% of an in-house baseline cohort, are categorized as germline polymorphisms.

3. **Reference-Track Artifacts** — False positives are mitigated by filtering out variants intersecting rCRS poly-N blocks, low-complexity repeats, or the standardized ENCODE hg38 genomic blacklist, consistent with core [mgatk](https://github.com/caleblareau/mgatk) recommendations.

---

## Somatic Nomination Rules

For each sample–variant tuple, scMOCHA partitions cellular records into cell-type annotations and evaluates the spatial distribution of the four evidence states. A variant is formally nominated as a somatic variant when it simultaneously satisfies three tunable operational boundaries:

1. **Recurrence** — Features ≥3 independent "colorful" cells within a minimum of one annotated cell type, with a global sample-level heteroplasmy between 5% and 95%.

2. **Coverage Breadth** — Requires sufficient reference alignment support ("black" cells) distributed across ≥6 or 7 distinct cell types. This ensures that the mutation's absence in other lineages represents true biological absence rather than localized data sparsity.

3. **Lineage Restriction** — Restricts wide-scale clonal distribution by verifying that the mutation is confined to few cell lineages:
   - *n*<sub>colorful</sub> ≤ 2 cell types for strictly lineage-restricted somatic events
   - *n*<sub>colorful</sub> ≤ 6 cell types for relaxed accumulation profiling

> Variants satisfying these distribution filters that inadvertently match known haplogroup paths or exceed population frequency thresholds are dynamically re-routed to a germline pool rather than reported as somatic events.

---

## Somatic-Mutation Hotspot

To profile cohort-wide mutational trends, the number of individuals carrying each distinct somatic variant is compiled and projected onto the 16,569-bp rCRS coordinate system using an interactive [gggenes](https://wilkox.org/gggenes/) lollipop layout stratified by major structural domain:

- **D-loop**
- **rRNA**
- **tRNA**
- **Protein-coding regions**
