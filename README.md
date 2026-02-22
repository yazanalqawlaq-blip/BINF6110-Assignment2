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

---

## Results

### Global Transcriptome Patterns

After filtering, 5,878 genes were retained for differential expression analysis. PCA of variance-stabilized counts revealed clear separation of samples by developmental stage along the first two principal components (**Figure 1**). PC1 captured the largest source of variance and separated Early biofilm samples from Thin and Mature samples, while PC2 distinguished Thin from Mature. Biological replicates clustered tightly within each stage, indicating high reproducibility and confirming that developmental stage is the dominant source of transcriptomic variation.

**Table 1.** Number of differentially expressed genes identified in each comparison (adjusted p < 0.05).

| Comparison | Upregulated | Downregulated | Total DE |
|-----------|-------------|---------------|----------|
| Thin vs. Early | 1,308 | 1,235 | 2,181* |
| Mature vs. Thin | 1,403 | 1,375 | 2,366* |
| Mature vs. Early | 1,737 | 1,625 | 2,958* |
| LRT (any stage change) | — | — | 3,712 |

*Totals include genes at all fold changes; up/down counts reflect adjusted p < 0.1 from DESeq2 summary with ashr shrinkage.

The number of DE genes increased with developmental distance, with the Mature vs. Early comparison showing the largest number of changes (2,958 genes, approximately 50% of all expressed genes). The LRT identified 3,712 genes (63% of expressed genes) with significant expression changes across any stage, confirming that velum maturation involves a widespread transcriptomic remodeling.

### Key Differentially Expressed Genes

**Table 2.** Top 20 most significantly differentially expressed genes in the Mature vs. Early comparison, ranked by adjusted p-value.

| Rank | ORF | Gene Name | log2FC | Adjusted p-value | Function |
|------|-----|-----------|--------|-------------------|----------|
| 1 | YJL052W | *TDH1* | −5.25 | 6.68 × 10⁻¹⁵¹ | Glyceraldehyde-3-phosphate dehydrogenase |
| 2 | YGL055W | *OLE1* | −4.63 | 5.97 × 10⁻¹³² | Δ9 fatty acid desaturase |
| 3 | YIR019C | *FLO11* | +5.53 | 1.83 × 10⁻¹³¹ | Cell-surface flocculin / adhesin |
| 4 | YGR087C | *PDC6* | −4.93 | 8.79 × 10⁻¹²⁴ | Pyruvate decarboxylase (minor isoform) |
| 5 | YHR094C | *HXT1* | −5.02 | 1.43 × 10⁻¹²¹ | Low-affinity hexose transporter |
| 6 | YNR073C | *MAN2* | +8.05 | 1.31 × 10⁻¹¹² | Mannosidase |
| 7 | YNR071C | — | +8.81 | 6.98 × 10⁻¹⁰⁰ | Uncharacterized ORF |
| 8 | YJR152W | *DAL5* | −3.41 | 4.15 × 10⁻⁹⁹ | Allantoate permease |
| 9 | YKL164C | *PIR1* | +3.94 | 1.54 × 10⁻⁹⁷ | Cell wall mannoprotein (PIR family) |
| 10 | YCR012W | *PGK1* | −4.16 | 1.12 × 10⁻⁸⁹ | Phosphoglycerate kinase |
| 11 | YNR072W | *HXT17* | +5.36 | 5.15 × 10⁻⁸⁷ | High-affinity hexose transporter |
| 12 | YOR273C | *TPO4* | −3.96 | 3.03 × 10⁻⁸³ | Polyamine transporter |
| 13 | YJL158C | *CIS3* | +4.40 | 7.18 × 10⁻⁷⁹ | Cell wall mannoprotein (PIR family) |
| 14 | YGR088W | *CTT1* | −4.43 | 2.72 × 10⁻⁷⁷ | Cytosolic catalase T |
| 15 | YPL106C | *SSE1* | +3.23 | 3.46 × 10⁻⁶⁸ | Hsp70 nucleotide exchange factor |
| 16 | YPL095C | *EEB1* | −3.77 | 1.40 × 10⁻⁶⁷ | Acyl-coenzymeA:ethanol O-acyltransferase |
| 17 | YGR060W | *ERG25* | −2.89 | 1.82 × 10⁻⁶⁵ | C-4 methyl sterol oxidase |
| 18 | YGR125W | *VSB1* | −2.31 | 3.50 × 10⁻⁶⁴ | Protein of unknown function |
| 19 | YJL212C | *OPT1* | −5.15 | 4.56 × 10⁻⁶¹ | Oligopeptide transporter |
| 20 | YOR247W | *SRL1* | +3.87 | 8.54 × 10⁻⁶⁰ | Cell wall mannoprotein |

### Functional Enrichment

GO Biological Process enrichment of DE genes (|log2FC| > 1, padj < 0.05) revealed that the most significantly enriched terms were transmembrane transport (188 genes, padj = 1.14 × 10⁻⁶), generation of precursor metabolites and energy (94 genes, padj = 3.00 × 10⁻⁵), and ATP metabolic process (31 genes, padj = 1.53 × 10⁻⁴) (**Figure 4**). The compareCluster analysis (**Figure 5**) showed that upregulated genes were enriched for mitochondrial translation and ribosome biogenesis, while downregulated genes were enriched for monocarboxylic acid metabolism, pyruvate metabolism, and organic acid metabolism.

GSEA analysis, which evaluates all genes ranked by fold change rather than just those passing a significance threshold, confirmed and extended the ORA findings. The most strongly enriched gene sets with negative normalized enrichment scores (NES) — indicating coordinated downregulation in mature biofilm — included purine compound catabolism (NES = −2.62), organophosphate catabolism (NES = −2.62), and monocarboxylic acid metabolism (NES = −2.40). Conversely, the most positively enriched gene sets included mitochondrial translation (NES = +2.38), cytoplasmic translation (NES = +1.89), and mitochondrion organization (NES = +2.14), all indicating upregulation in mature biofilm.

---

<img width="2100" height="1800" alt="Image" src="https://github.com/user-attachments/assets/76fc0b3c-0df0-4890-bf85-abc4207fccc1" />

**Figure 1. PCA of Yeast Biofilm Samples by Development Stage.** Principal component analysis of variance-stabilized expression data for 9 samples across three stages of velum development. PC1 and PC2 capture the majority of expression variance and cleanly separate all three developmental stages, with biological replicates clustering tightly within each stage.


<img width="2400" height="2100" alt="Image" src="https://github.com/user-attachments/assets/23853b51-e9c3-4130-b933-9f43b3f6e238" />

**Figure 2. Volcano Plot of Differentially Expressed Genes (Mature vs. Early Biofilm).** Each point represents a gene, plotted by log2 fold change (x-axis) and statistical significance (y-axis). Red points indicate significantly upregulated genes in mature biofilm (log2FC > 1, padj < 0.05), blue points indicate downregulated genes, and gray points are not significant. Dashed lines indicate the fold change (±1) and significance (padj = 0.05) thresholds. The top 15 genes by adjusted p-value are labeled. *FLO11*, the key biofilm adhesin, is among the most significantly upregulated genes, while glycolytic enzymes (*TDH1*, *PGK1*, *PDC6*) and the low-affinity hexose transporter *HXT1* are among the most strongly downregulated.


<img width="2400" height="3000" alt="Image" src="https://github.com/user-attachments/assets/b4a89320-2700-4d69-9fd9-9c5bc2eaa61f" />

**Figure 3. Heatmap of the Top 30 Differentially Expressed Genes Across Biofilm Stages.** Row-scaled expression (VST) of the 30 most significant genes identified by the likelihood ratio test. Hierarchical clustering of columns separates samples cleanly by developmental stage. Two major gene clusters are visible: one comprising genes highly expressed in Early biofilm (including glycolytic genes *TDH1*, *PGK1*, and the fatty acid desaturase *OLE1*) and another comprising genes upregulated in Mature biofilm (including *FLO11*, cell wall mannoproteins *PIR1*, *CIS3*, and *SRL1*, and the hexose transporter *HXT17*).


<img width="2700" height="2100" alt="Image" src="https://github.com/user-attachments/assets/f2be467d-fc1f-47f1-94ac-a2eaa6509f10" />

**Figure 4. GO Biological Process Enrichment of DE Genes (Mature vs. Early).** Dotplot showing the top 15 enriched GO Biological Process terms among genes differentially expressed between mature and early biofilm. Point size reflects gene count and color reflects adjusted p-value. Transmembrane transport, energy metabolism, and organic acid metabolism are the most prominent enriched processes.


<img width="3000" height="2400" alt="Image" src="https://github.com/user-attachments/assets/87fcad5e-44f8-4f75-8439-8e995e8da48d" />

**Figure 5. GO Biological Process Enrichment Comparing Upregulated vs. Downregulated Genes.** compareCluster dotplot showing the top 10 enriched GO BP terms separately for genes upregulated and downregulated in mature versus early biofilm. This directional analysis reveals that mitochondrial translation and ribosome assembly processes are enriched among upregulated genes, while metabolic and catabolic processes dominate among downregulated genes.


<img width="2700" height="2100" alt="Image" src="https://github.com/user-attachments/assets/15471344-0db0-4758-834a-3c26303ffe4d" />

**Figure 6. KEGG Pathway Enrichment of DE Genes.** Dotplot showing enriched KEGG pathways among significantly differentially expressed genes.


<img width="2400" height="1800" alt="Image" src="https://github.com/user-attachments/assets/180d57e0-0574-46d7-9500-9a93a840c433" />

**Figure 7. MA Plot (Mature vs. Early Biofilm).** Mean expression (x-axis) versus log2 fold change (y-axis) for all genes. Red points indicate significantly DE genes (padj < 0.05). The symmetric spread of points above and below zero indicates balanced upregulation and downregulation. Highly expressed genes (right side) show more moderate fold changes, consistent with the biological expectation that housekeeping genes are less likely to be dramatically altered.

---

## Discussion

### FLO11 and the Structural Basis of Velum Formation

The most biologically significant finding in this analysis is the strong upregulation of *FLO11* (log2FC = +5.53, padj = 1.83 × 10⁻¹³¹), making it the third most statistically significant gene in the entire dataset. *FLO11* encodes the primary cell-surface adhesin required for velum formation in flor yeasts (Zara et al., 2005). The Flo11 protein is a GPI-anchored glycoprotein that mediates both cell–cell adhesion and cell–surface hydrophobicity, both of which are essential for cells to aggregate at the air–liquid interface (Dranginis et al., 2007). The 46-fold increase in *FLO11* expression between Early and Mature biofilm is consistent with the findings of Mardanov et al. (2020), who reported that *FLO11* was the only flocculin gene strongly upregulated during velum maturation, while other flocculins (*FLO1*, *FLO5*, *FLO9*, *FLO10*) remained unchanged or were moderately downregulated.

This specific upregulation of *FLO11* — but not other flocculins — is biologically meaningful because it reflects the unique structural requirements of velum formation. While Flo1p and other flocculins primarily drive cell–cell aggregation in liquid (flocculation), Flo11p uniquely confers the cell-surface hydrophobicity needed to maintain a stable buoyant biofilm at the air–liquid interface (Váchová et al., 2021). This is consistent with the observation that flor strains carry specific promoter modifications in the *ICR1* noncoding RNA region upstream of *FLO11* that enhance its transcription, suggesting strong selective pressure for high *FLO11* expression in the velum niche (Fidalgo et al., 2006; Eldarov et al., 2018).

Alongside *FLO11*, several cell wall mannoproteins were strongly upregulated: *PIR1* (log2FC = +3.94), *CIS3*/PIR4 (log2FC = +4.40), and *SRL1* (log2FC = +3.87). The PIR (Proteins with Internal Repeats) family contributes to cell wall integrity and cross-linking, and their co-upregulation with *FLO11* suggests a coordinated remodeling of the cell wall surface during biofilm maturation (Essen et al., 2020). This strengthened cell wall is likely important for maintaining structural integrity of the velum, which must support itself as a floating mat while being exposed to mechanical stress at the air–liquid interface.

### Metabolic Reprogramming: From Fermentation to Oxidative Metabolism

The most striking pattern among downregulated genes is the suppression of glycolytic enzymes and fermentation genes. *TDH1* (glyceraldehyde-3-phosphate dehydrogenase, log2FC = −5.25), *PGK1* (phosphoglycerate kinase, log2FC = −4.16), and *PDC6* (pyruvate decarboxylase, log2FC = −4.93) are all core components of the glycolytic pathway and ethanol fermentation. Their strong downregulation in mature biofilm reflects the well-documented metabolic shift that occurs during velum development: as fermentable sugars are depleted and cells gain access to oxygen at the biofilm surface, they transition from fermentative to oxidative metabolism, using ethanol and acetaldehyde as primary carbon sources (Alexandre, 2013; Moreno-García et al., 2018).

This metabolic interpretation is further supported by the hexose transporter expression patterns. *HXT1*, a low-affinity glucose transporter active when glucose is abundant, was one of the most strongly downregulated genes (log2FC = −5.02). In contrast, *HXT17*, a high-affinity hexose transporter expressed under glucose-limiting conditions, was strongly upregulated (log2FC = +5.36). This switch from low-affinity to high-affinity sugar transport is a classic marker of the glucose repression/derepression transition in yeast and confirms that mature biofilm cells are operating in a carbon-limited, non-fermentative metabolic state (Özcan and Johnston, 1999).

The GSEA results reinforce this picture at the pathway level: monocarboxylic acid metabolism, organophosphate catabolism, and lipid biosynthesis were all coordinately downregulated (negative NES values), while mitochondrial translation (NES = +2.38) and mitochondrion organization (NES = +2.14) were upregulated. The activation of mitochondrial gene expression is entirely consistent with the shift to respiratory metabolism, as cells require increased mitochondrial capacity to oxidize ethanol through the TCA cycle and electron transport chain rather than producing it through fermentation (Moreno-García et al., 2018).

### Additional Genes of Biological Interest

Several other highly significant DE genes provide additional biological insight into the velum maturation process. *OLE1* (Δ9 fatty acid desaturase, log2FC = −4.63) catalyzes the introduction of a double bond into saturated fatty acids and is essential for membrane fluidity. Its strong downregulation may reflect changes in membrane lipid composition as cells adapt to the oxygen-rich environment at the biofilm surface, where membrane remodeling is part of the stress response (Alexandre, 2013). *EEB1* (ethanol acyltransferase, log2FC = −3.77) contributes to the synthesis of ethyl esters, which are important flavor compounds in wine; its downregulation at the mature stage is consistent with the shift away from fermentation-associated volatile production. *CTT1* (cytosolic catalase T, log2FC = −4.43) protects against oxidative stress; its unexpected downregulation may indicate that mature biofilm cells rely on alternative antioxidant systems or have reduced intracellular reactive oxygen species due to metabolic changes at the biofilm surface.

*MAN2* (mannosidase, log2FC = +8.05), the gene with the largest positive fold change among annotated genes, is involved in glycoprotein processing and cell wall mannan metabolism. Its dramatic upregulation, together with the PIR family mannoproteins, further underscores that cell wall restructuring is one of the most prominent transcriptomic signatures of velum maturation. The uncharacterized ORF *YNR071C* (log2FC = +8.81) is located adjacent to *HXT17* and *MAN2* on chromosome XIV, and its very high fold change suggests it may play a currently unknown role in the biofilm developmental program; this is consistent with Mardanov et al. (2020), who noted that several highly expressed uncharacterized genes in mature biofilm showed stress-associated expression patterns.

### Comparison with Published Results

The overall patterns I observed are highly consistent with Mardanov et al. (2020), despite methodological differences (Salmon/DESeq2 vs. HISAT2/DESeq). Both analyses identified *FLO11* as the most prominently upregulated adhesin, confirmed the downregulation of glycolytic genes, and revealed the shift in hexose transporter expression. The GO enrichment results also align closely: transmembrane transport, energy metabolism, and organic acid metabolism were the dominant enriched categories in both analyses. The consistency of results across different quantification and analysis pipelines strengthens confidence that the observed expression changes reflect genuine biology rather than method-specific artifacts. One difference is that my analysis detected a larger total number of DE genes (2,958 at padj < 0.05 for Mature vs. Early), which may reflect the improved sensitivity of DESeq2 with ashr shrinkage compared to the older DESeq methodology used in the original study.

### Limitations

This analysis has several limitations worth noting. The dataset includes only three biological replicates per stage, which limits statistical power for detecting subtle expression changes and may underestimate the true number of DE genes for comparisons between adjacent stages (Thin vs. Early, Mature vs. Thin). The use of single-end reads, while sufficient for gene-level quantification, prevents analysis of alternative splicing events that might accompany the developmental transition. Additionally, the three time points (38, 83, 109 days) capture snapshots of a continuous process, and the large temporal gaps between sampling points mean that the timing and dynamics of individual gene expression changes cannot be precisely resolved. Finally, KEGG pathway analysis initially failed due to identifier mapping issues between ORF names and KEGG gene IDs, highlighting a practical challenge when working with organism-specific annotation databases.

## Conclusion

This reanalysis of the Mardanov et al. (2020) yeast biofilm RNA-seq dataset confirms that velum maturation involves a profound transcriptomic reprogramming affecting the majority of expressed genes. The central findings — strong upregulation of *FLO11* and cell wall mannoproteins, downregulation of glycolytic and fermentative genes, and a switch from low-affinity to high-affinity sugar transport — together paint a coherent picture of yeast cells transitioning from a fermentative lifestyle to a surface-attached, oxygen-utilizing biofilm state. The functional enrichment analyses further support this interpretation, revealing coordinated activation of mitochondrial and translational machinery alongside suppression of fermentative metabolic pathways. These results demonstrate that modern lightweight quantification pipelines (Salmon + DESeq2 + clusterProfiler) can effectively recapitulate and extend the findings of the original study, while providing additional biological insight through comprehensive functional annotation.

---

## References

Alexandre, H. (2013). Flor yeasts of *Saccharomyces cerevisiae* — their ecology, genetics and metabolism. *International Journal of Food Microbiology*, 167(2), 269–275. https://doi.org/10.1016/j.ijfoodmicro.2013.08.021

Anders, S., & Huber, W. (2010). Differential expression analysis for sequence count data. *Genome Biology*, 11(10), R106. https://doi.org/10.1186/gb-2010-11-10-r106

Brückner, S., & Mösch, H.-U. (2012). Choosing the right lifestyle: adhesion and development in *Saccharomyces cerevisiae*. *FEMS Microbiology Reviews*, 36(1), 25–58. https://doi.org/10.1111/j.1574-6976.2011.00275.x

Dranginis, A. M., Rauceo, J. M., Coronado, J. E., & Lipke, P. N. (2007). A biochemical guide to yeast adhesins: glycoproteins for social and antisocial occasions. *Microbiology and Molecular Biology Reviews*, 71(2), 282–294. https://doi.org/10.1128/MMBR.00037-06

Eldarov, M. A., Beletsky, A. V., Tanashchuk, T. N., Kishkovskaya, S. A., Ravin, N. V., & Mardanov, A. V. (2018). Whole-genome analysis of three yeast strains used for production of sherry-like wines revealed genetic traits specific to flor yeasts. *Frontiers in Microbiology*, 9, 965. https://doi.org/10.3389/fmicb.2018.00965

Essen, L.-O., Vogt, M. S., & Mösch, H.-U. (2020). Diversity of GPI-anchored fungal adhesins. *Biological Chemistry*, 401(12), 1389–1405. https://doi.org/10.1515/hsz-2020-0199

Fidalgo, M., Barrales, R. R., Ibeas, J. I., & Jimenez, J. (2006). Adaptive evolution by mutations in the *FLO11* gene. *Proceedings of the National Academy of Sciences*, 103(30), 11228–11233. https://doi.org/10.1073/pnas.0601713103

Legras, J.-L., Erny, C., & Charpentier, C. (2016). Population structure and comparative genome hybridization of European flor yeast reveal a unique group of *Saccharomyces cerevisiae* strains with few gene duplications in their genome. *PLOS ONE*, 9(10), e108089. https://doi.org/10.1371/journal.pone.0108089

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8

Mardanov, A. V., Eldarov, M. A., Beletsky, A. V., Tanashchuk, T. N., Kishkovskaya, S. A., & Ravin, N. V. (2020). Transcriptome profile of yeast strain used for biological wine aging revealed dynamic changes of gene expression in course of flor development. *Frontiers in Microbiology*, 11, 538. https://doi.org/10.3389/fmicb.2020.00538

Moreno-García, J., Mauricio, J. C., Moreno, J., & García-Martínez, T. (2018). Differential proteome analysis of a flor yeast strain under biofilm formation. *International Journal of Molecular Sciences*, 18(3), 720. https://doi.org/10.3390/ijms18030720

Özcan, S., & Johnston, M. (1999). Function and regulation of yeast hexose transporters. *Microbiology and Molecular Biology Reviews*, 63(3), 554–569. https://doi.org/10.1128/MMBR.63.3.554-569.1999

Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods*, 14(4), 417–419. https://doi.org/10.1038/nmeth.4197

Purevdorj-Gage, B., Orr, M. E., Stoodley, P., Sheehan, K. B., & Hyman, L. E. (2007). The role of FLO11 in *Saccharomyces cerevisiae* biofilm development in a laboratory based flow-cell system. *FEMS Yeast Research*, 7(3), 372–379. https://doi.org/10.1111/j.1567-1364.2006.00189.x

Soneson, C., Love, M. I., & Robinson, M. D. (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research*, 4, 1521. https://doi.org/10.12688/f1000research.7563.2

Stephens, M. (2016). False discovery rates: a new deal. *Biostatistics*, 18(2), 275–294. https://doi.org/10.1093/biostatistics/kxw041

Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., ... & Mesirov, J. P. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proceedings of the National Academy of Sciences*, 102(43), 15545–15550. https://doi.org/10.1073/pnas.0506580102

Váchová, L., Šťovíček, V., Hlaváček, O., Chernyavskiy, O., Štěpánek, L., Kubínová, L., & Palková, Z. (2021). *FLO11*, a developmental gene conferring impressive adaptive plasticity to the yeast *Saccharomyces cerevisiae*. *Pathogens*, 10(11), 1509. https://doi.org/10.3390/pathogens10111509

Wu, T., Hu, E., Xu, S., Chen, M., Guo, P., Dai, Z., Feng, T., Zhou, L., Tang, W., Zhan, L., Fu, X., Liu, S., Bo, X., & Yu, G. (2021). clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. *The Innovation*, 2(3), 100141. https://doi.org/10.1016/j.xinn.2021.100141

Zara, S., Bakalinsky, A. T., Zara, G., Pirino, G., Demontis, M. A., & Budroni, M. (2005). *FLO11*-based model for air-liquid interfacial biofilm formation by *Saccharomyces cerevisiae*. *Applied and Environmental Microbiology*, 71(6), 2934–2939. https://doi.org/10.1128/AEM.71.6.2934-2939.2005
