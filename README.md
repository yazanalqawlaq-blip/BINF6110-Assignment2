#  Transcriptome Analysis of Yeast Biofilm (Velum) Development During Biological Wine Aging

**Course:** BINF*6110 - Bioinformatics

**Author:** Yazan Alqawlaq

**Date:** March 2026

**Repository:** https://github.com/yazanalqawlaq-blip/BINF6110-Assignment2

---

## Introduction

Flor strains of *Saccharomyces cerevisiae* are specialized yeasts responsible for the biological aging of sherry-type wines, a process in which yeast cells form a floating biofilm called a velum (or flor) at the air–liquid interface of fortified wine (Alexandre, 2013). This velum is a major adaptive phenotype that allows flor yeasts to survive in the harsh winemaking environment, which is characterized by high ethanol concentrations (10–15% v/v), nutrient depletion, and low oxygen availability in the liquid phase (Legras et al., 2016). By forming a biofilm at the surface, yeast cells gain direct access to atmospheric oxygen, enabling a metabolic shift from fermentation to oxidative metabolism of ethanol and other non-fermentable carbon sources (Moreno-García et al., 2018). This metabolic transition drives the characteristic biochemical and sensory changes associated with biological wine aging, including the accumulation of acetaldehyde, acetals, and other volatile compounds that define the flavor profile of sherry wines (Mardanov et al., 2020).

Velum formation depends critically on the cell-surface adhesin encoded by *FLO11*, a member of the flocculin gene family in *S. cerevisiae* (Dranginis et al., 2007). Unlike other flocculins (Flo1p, Flo5p, Flo9p, Flo10p), which primarily mediate calcium-dependent cell–cell aggregation (flocculation), Flo11p is unique in enabling a broad range of adhesion-dependent phenotypes, including biofilm formation, velum formation, invasive growth, and pseudohyphal development (Zara et al., 2005; Purevdorj-Gage et al., 2007). The regulation of *FLO11* expression is exceptionally complex, involving one of the largest promoters in the yeast genome and integrating signals from the MAPK, cAMP-PKA, and TOR signaling pathways as well as epigenetic mechanisms involving noncoding RNAs (Brückner and Mösch, 2012; Váchová et al., 2021). Flor yeast strains often carry specific promoter deletions in the *ICR1* noncoding RNA region that constitutively enhance *FLO11* transcription, which is thought to be a key genetic adaptation underlying their capacity for velum formation (Fidalgo et al., 2006; Eldarov et al., 2018).

Mardanov et al. (2020) generated the first comprehensive RNA-seq transcriptome dataset for a flor yeast strain across multiple stages of velum development under experimental winemaking conditions. Using the industrial flor strain I-329, they collected samples at three developmental stages spanning 71 days: early biofilm (38 days post-inoculation), thin biofilm (83 days), and mature biofilm (109 days), with three biological replicates per stage. Their analysis revealed dynamic changes in gene expression involving carbon metabolism, cell wall biogenesis, stress response, and adhesion genes, but their differential expression analysis relied on the older DESeq methodology (Anders and Huber, 2010) with a reference-guided alignment approach. In this assignment, I reanalyze the Mardanov et al. (2020) RNA-seq dataset using a modern lightweight quantification pipeline based on Salmon and DESeq2, followed by functional enrichment analysis with clusterProfiler. This approach allows me to independently evaluate the transcriptomic changes occurring during velum maturation and to examine how specific genes and pathways contribute to the biological transition from early to mature biofilm.

I chose Salmon (Patro et al., 2017) for transcript quantification because it uses a quasi-mapping algorithm that avoids the computational cost of full read alignment while providing accurate abundance estimates with built-in correction for GC-content bias and sequence-specific bias. Compared to traditional alignment-based approaches such as HISAT2/featureCounts or STAR/HTSeq, Salmon is substantially faster and has been shown to yield more accurate transcript-level estimates, particularly for genes belonging to families with similar sequences (Patro et al., 2017). One limitation of alignment-free methods is that they do not produce a BAM file for visual inspection of read pileups, which means some types of quality assessment (such as checking for unexpected splice patterns or contamination at specific loci) are not possible without a separate alignment step. However, for standard differential expression workflows, the speed and accuracy advantages of Salmon are well established (Soneson et al., 2015). I used DESeq2 (Love et al., 2014) for differential expression because it implements a generalized linear model framework with shrinkage estimation of dispersions, which is designed for RNA-seq experiments with small sample sizes and provides robust statistical inference. Compared to edgeR, which uses a different approach to dispersion estimation and can be more liberal in calling DE genes, DESeq2 tends to be more conservative and produces well-calibrated false discovery rates (Love et al., 2014). I used clusterProfiler (Wu et al., 2021) for functional enrichment because it provides a unified framework for Gene Ontology and KEGG over-representation analysis as well as gene set enrichment analysis (GSEA), and it integrates well with organism-specific annotation databases for yeast.

---

### Project Structure
```
BINF6110-Assignment2/
├── input_data/              # Raw reads and reference files (gitignored)
├── output_files/
│   ├── fastqc/              # FastQC quality reports
│   ├── salmon_quant/        # Salmon quantification output
│   └── R_analysis/          # DESeq2, enrichment, and figures
│       ├── DE_*.csv          # Differential expression results
│       ├── GO_*_enrichment.csv
│       ├── GSEA_GO_BP.csv
│       └── figures/          # Publication-quality figures
├── scripts/                 # Full reproducible pipeline (01–08)
└── README.md                # Complete documentation
```

---

## Methods

### Data Acquisition

RNA-seq data from Mardanov et al. (2020) was obtained from the NCBI Sequence Read Archive (BioProject PRJNA592304). Nine single-end Illumina RNA-seq libraries corresponding to three biological replicates at each of three velum development stages were downloaded: early biofilm at 38 days (SRR10551665, SRR10551664, SRR10551663), thin biofilm at 83 days (SRR10551662, SRR10551661, SRR10551660), and mature biofilm at 109 days (SRR10551659, SRR10551658, SRR10551657). Data was downloaded using `prefetch` and converted to FASTQ format using `fasterq-dump` (SRA Toolkit v3.2.1) with `--split-3`, followed by compression with `gzip`. All downloads were performed on the Narval login node as compute nodes lack internet access.

The *S. cerevisiae* R64 reference genome (accession GCF_000146045.2), transcriptome, and gene annotation (GTF) were downloaded from NCBI RefSeq.

### Quality Control

Read quality was assessed using FastQC v0.12.1 on all nine samples. All samples showed high quality scores (median Phred > 30 across all positions) with no adapter contamination, so no trimming was performed.

### Transcript Quantification

Transcript quantification was performed using Salmon v1.10.1 (Patro et al., 2017) in quasi-mapping mode. A decoy-aware index was constructed following the recommended approach from BINF*6110 Lecture 6, where the full genome sequence was used as a decoy to reduce spurious mapping of reads derived from intergenic or intronic regions. The genome FASTA was concatenated after the transcriptome FASTA to produce a combined `gentrome.fna` file, and chromosome names were extracted to generate a `decoys.txt` file. The Salmon index was built with default k-mer size (k=31).

Quantification was performed on all nine samples with the following parameters: automatic library type detection (`-l A`), selective alignment validation (`--validateMappings`), GC-content bias correction (`--gcBias`), and estimated fragment length parameters for single-end reads (`--fldMean 180 --fldSD 20`). The `--fldMean` and `--fldSD` flags are required for single-end libraries because Salmon cannot infer the fragment length distribution empirically without paired-end information; the values of 180 bp mean and 20 bp standard deviation represent reasonable estimates for standard Illumina RNA-seq library preparations. Each sample yielded 4.4–7.5 million mapped reads.

### Differential Expression Analysis

Salmon quantification files were imported into R (v4.5.1) using tximport v1.36.0 (Soneson et al., 2015), which aggregated transcript-level estimates to gene-level counts using a transcript-to-gene mapping derived from the GTF annotation via GenomicFeatures. Tximport is preferred over naive count summation because it accounts for transcript length differences and the uncertainty of multi-mapped reads when collapsing to gene level (Soneson et al., 2015).

Differential expression analysis was performed with DESeq2 v1.48.1 (Love et al., 2014). A count matrix of 6,046 genes across 9 samples was constructed, and genes with fewer than 10 total reads across all samples were filtered, leaving 5,878 genes for analysis. The DESeq2 model was fit with the design formula `~ stage`, where stage is a three-level factor (Early, Thin, Mature).

Three pairwise comparisons were performed using Wald tests: Thin vs. Early, Mature vs. Thin, and Mature vs. Early. Log2 fold change estimates were shrunk using the adaptive shrinkage estimator `ashr` (Stephens, 2016) to obtain more accurate effect sizes, particularly for genes with low counts or high dispersion. In addition, a likelihood ratio test (LRT) was performed comparing the full model (`~ stage`) to the reduced model (`~ 1`) to identify genes whose expression changed significantly across any stage of biofilm development, regardless of the direction or timing of the change. Genes were considered significantly differentially expressed at an adjusted p-value threshold of 0.05 (Benjamini–Hochberg correction).

### Functional Enrichment Analysis

Functional enrichment was performed using clusterProfiler v4.16.0 (Wu et al., 2021) with the *S. cerevisiae* annotation package org.Sc.sgd.db. Gene identifiers (yeast ORF names, e.g. YIR019C) were mapped to standard gene names and Entrez IDs using `bitr()`.

Over-representation analysis (ORA) was performed on significantly differentially expressed genes from the Mature vs. Early comparison (|log2FC| > 1, adjusted p < 0.05) using `enrichGO()` for Gene Ontology Biological Process (BP), Molecular Function (MF), and Cellular Component (CC) terms, with the full set of 5,878 expressed genes as the background universe. A `compareCluster()` analysis was also performed to separately examine GO BP enrichment in upregulated (905 genes) versus downregulated (745 genes) gene sets, which provides biological insight into the directionality of pathway changes.

Gene set enrichment analysis (GSEA) was performed using `gseGO()` with all genes ranked by log2 fold change from the Mature vs. Early comparison. Unlike ORA, which requires an arbitrary significance cutoff, GSEA evaluates whether predefined gene sets are enriched at the extremes of a ranked list, making it a complementary approach that can detect subtle but coordinated expression shifts (Subramanian et al., 2005).

All enrichment results were filtered at adjusted p-value < 0.05 (BH correction) and q-value < 0.2.

### Visualization

Seven figures were generated to comprehensively characterize the transcriptomic landscape of velum development. Principal component analysis (PCA) was performed on variance-stabilized counts (VST) to assess global sample clustering (Figure 1). A volcano plot was generated for the Mature vs. Early comparison with the top 15 most significant genes labeled by common gene name (Figure 2). A heatmap of the top 30 genes from the LRT analysis was generated using pheatmap with row-scaled expression values and hierarchical clustering of both genes and samples (Figure 3). GO Biological Process enrichment was visualized as a dotplot showing the top 15 enriched terms (Figure 4). A compareCluster dotplot contrasted enriched GO terms between upregulated and downregulated genes (Figure 5). KEGG pathway enrichment was displayed as a dotplot (Figure 6). An MA plot was generated to show the relationship between mean expression and log2 fold change (Figure 7).

### Computational Environment

Shell-based analyses (data download, FastQC, Salmon indexing and quantification) were executed on the Compute Canada Narval HPC cluster using SLURM job scheduling. Bioinformatics tools were loaded as environment modules (`StdEnv/2023`, `sra-toolkit/3.2.1`, `fastqc/0.12.1`, `salmon/1.10.1`). R-based analyses (DESeq2, clusterProfiler, visualization) were performed locally in RStudio. All scripts, parameters, and outputs are documented in the project GitHub repository.

### Software Versions

| Tool | Version | Purpose |
|------|---------|---------|
| SRA Toolkit | 3.2.1 | Data download and FASTQ conversion |
| FastQC | 0.12.1 | Read quality assessment |
| Salmon | 1.10.1 | Transcript quantification |
| R | 4.5.1 | Statistical computing environment |
| tximport | 1.36.0 | Import Salmon counts to gene level |
| DESeq2 | 1.48.1 | Differential expression analysis |
| ashr | 2.2-63 | Log2 fold change shrinkage |
| clusterProfiler | 4.16.0 | GO and KEGG enrichment analysis |
| org.Sc.sgd.db | 3.20.0 | *S. cerevisiae* gene annotation |
| pheatmap | 1.0.13 | Heatmap visualization |
| ggplot2 | 3.5.2 | Figure generation |

