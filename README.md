# Population Genomics Pipeline for the Critically Endangered Spoon-billed Sandpiper (*Calidris pygmaea*)

This repository provides a fully automated, reproducible computational workflow implemented in **Snakemake** for the comprehensive analysis of whole-genome sequencing (WGS) data from the spoon-billed sandpiper (*Calidris pygmaea*), a Critically Endangered migratory shorebird. The pipeline encompasses the entire analytical chain — from raw variant calling output to demographic inference, estimation of genetic load, inbreeding assessment, and forward-in-time evolutionary simulations.

## Workflow Overview and Analytical Modules

The pipeline is structured into eight logically distinct modules, ensuring modularity, computational efficiency, and full reproducibility.

### 1. Sex Chromosome Identification

Given the scaffold-level assembly of the reference genome, Z and W sex chromosomes are not reliably annotated.  
- **Approach**: Homology-based identification using reference Z (zebra finch, *Taeniopygia guttata*) and W (chicken, *Gallus gallus*) chromosomes retrieved from NCBI.  
- **Tool**: `BLASTn`  
- **Output**: Ranked list of scaffolds showing significant homology to sex chromosomes, enabling their selective inclusion or exclusion in downstream analyses.

### 2. Variant Quality Control and Filtering

Transformation of raw multi-sample VCF files into analysis-ready datasets.  
- **Tools**: `bcftools`, `vcftools`  
- **Filtering strategy**: Removal of indels, multi-allelic sites, and low-quality variants (based on QUAL and other standard metrics).  
- **Dual-track design**:  
  - **Full VCF** — retains all high-quality SNPs; used for mutational load and nucleotide diversity (π) estimation.  
  - **Pruned VCF** — linkage disequilibrium (LD) pruned dataset; used for PCA, admixture analysis, and demographic modeling in dadi.

### 3. Population Structure Analysis

Assessment of genetic homogeneity across the sampled individuals (n = 22).  
- **Principal Component Analysis (PCA)**: dimensionality reduction and visualization performed with `PLINK`.  
- **Sparse Non-negative Matrix Factorization (sNMF)**: ancestry proportion estimation and cross-entropy evaluation to infer optimal number of ancestral populations (K) using the LEA package in R.  
- **Technical solution**: On-the-fly scaffold-to-chromosome renaming (mapping to chromosome 1) to circumvent PLINK 1.9 limitations with large numbers of unplaced scaffolds.

### 4. Demographic History Reconstruction (dadi)

Inference of historical effective population size changes.  
- **Method**: Allele frequency spectrum (SFS)-based modeling.  
- **Focal model**: Bottleneck followed by recovery.  
- **Implementation**: automated parameter optimization (θ = 4Nₑμ, timing of events, recovery rate); best-fit parameters are exported for downstream forward simulations.

### 5. Estimation of Mutational Load

Quantification of deleterious variant accumulation.  
- **Tool**: `SnpEff` (functional annotation using provided GFF)  
- **Variant classification**:  
  - Synonymous (neutral benchmark)  
  - Missense (mildly deleterious)  
  - Loss-of-function (LoF; strongly deleterious)  
- **Summary metric**: Ratio of deleterious (missense + LoF) to synonymous variants per individual.

### 6. Forward-in-Time Evolutionary Simulation (SLiM 4)

Testing the hypothesis of genetic purging versus deleterious mutation accumulation following a severe bottleneck.  
- **Approach**: Individual-based forward simulations.  
- **Integration**: Realistic demographic parameters (Nₑ trajectory, timing) derived directly from the best-fit dadi model.  
- **Aim**: Evaluate long-term fate of deleterious alleles under the inferred demographic history.

### 7. Runs of Homozygosity (ROH) and Inbreeding

Detection and characterization of autozygous segments.  
- **Tool**: `PLINK --homozyg`  
- **Metrics**:  
  - Inbreeding coefficient *F*<sub>ROH</sub>  
  - Length-based classification (short → ancient inbreeding; long → recent mating among relatives)  
- **Visualization**: Distribution of ROH lengths and per-individual summaries.

### 8. Genomic Window-based Scans

Spatial heterogeneity of genetic diversity across the genome.  
- **Metrics**: Nucleotide diversity (π) and Tajima's D  
- **Resolution**: Sliding windows of 50 kb  
- **Purpose**: Identification of putative selective sweep regions or loci associated with key adaptations (e.g., bill morphology).

## Technology Stack

- **Workflow engine**: Snakemake  
- **Scripting languages**: Python 3, R, Bash  
- **Core bioinformatics tools**: bcftools, vcftools, PLINK 1.9, BLAST+, SnpEff, SLiM 4, dadi  
- **Visualization libraries**: Matplotlib, Seaborn, LEA (R)
