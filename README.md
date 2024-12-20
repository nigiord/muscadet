> [!NOTE]
> This R package is currently in development and therefore is not functional yet

# muscadet <img src="man/figures/logo.png" align="right" height="139" alt="" />

### *multiomics single-cell copy number alterations detection*

**R package for the identification of copy number alterations (CNAs) in cancer cells from single-cell multiomics data.**


> Marie Denoulet<sup>1</sup>, Mia Cherkaoui<sup>1</sup>, Nils Giordano<sup>1</sup>, Robin Lanée<sup>1</sup>, Elise Douillard<sup>1,2</sup>, Magali Devic<sup>1,2</sup>, Florence Magrangeas<sup>1,2</sup>, Stéphane Minvielle<sup>1,2</sup>, Céline Vallot<sup>3,4</sup>, Eric Letouzé<sup>1,2</sup>

> <sup>1</sup>Nantes Université, INSERM, CNRS, Université d'Angers, CRCI2NA, Nantes, France. <sup>2</sup>University Hospital Hôtel-Dieu, Nantes, France. <sup>3</sup>CNRS UMR3244, Institut Curie, PSL University, Paris, France. <sup>4</sup>Translational Research Department, Institut Curie, PSL University, Paris, France

### Install the latest version directly from GitHub
```r
library(devtools)
install_github("ICAGEN/muscadet")
```

### Use muscadet with example data
```r
library(muscadet)

# Create muscomic objects
atac <- CreateMuscomicObject(
    type = "ATAC",
    mat_counts = mat_counts_atac_tumor,
    allele_counts = allele_counts_atac_tumor,
    features = peaks)
rna <- CreateMuscomicObject(
    type = "RNA",
    mat_counts = mat_counts_rna_tumor,
    allele_counts = allele_counts_rna_tumor,
    features = genes)
atac_ref <- CreateMuscomicObject(
    type = "ATAC",
    mat_counts = mat_counts_atac_ref,
    allele_counts = allele_counts_atac_ref,
    features = peaks)
rna_ref <- CreateMuscomicObject(
    type = "RNA",
    mat_counts = mat_counts_rna_ref,
    allele_counts = allele_counts_rna_ref,
    features = genes)

# Create raw muscadet objects
muscadet <- CreateMuscadetObject(
    omics = list(ATAC = atac, RNA = rna),
    bulk.lrr = bulk_lrr,
    bulk.label = "WGS",
    genome = "hg38")
muscadet_ref <- CreateMuscadetObject(
    omics = list(ATAC = atac_ref, RNA = rna_ref),
    genome = "hg38")

# Compute log R ratios from scATAC-seq read counts
muscadet <- computeLogRatio(
    x = muscadet,
    reference = muscadet_ref,
    omic = "ATAC",
    method = "ATAC",
    minReads = 1,
    minPeaks = 1)

# Compute log R ratios from scRNA-seq read counts
muscadet <- computeLogRatio(
    x = muscadet,
    reference = muscadet_ref,
    omic = "RNA",
    method = "RNA",
    refReads = 2)

# Cluster cells based on log R ratio data
muscadet <- clusterMuscadet(
    muscadet,
    dist_method = "euclidean",
    hclust_method = "ward.D",
    k_range = 2:5)

# You can also load the example muscadet object
load(muscadet_obj)

# Plot heatmap of clustering
heatmapMuscadet(
    muscadet,
    k = 3,
    show_missing = FALSE,
    filename = "heatmap_muscadet_k3_commoncells.png",
    title = "Example sample (k=3)")
heatmapMuscadet(
    muscadet,
    k = 3,
    filename = "heatmap_muscadet_k3_allcells.png",
    title = "Example sample (k=3)")

# Plot Silhouette widths for cluster validation
plotSil(muscadet, k = 3)

# Plot cluster validation indexes
plotIndexes(muscadet)
plotIndexes(muscadet, index = "silhouette")

# Assign the chosen cluster assignment
muscadet <- assignClusters(muscadet, k = 3)

# Merge all counts from all omics from both sample and reference
muscadet <- mergeCounts(muscadet, muscadet_ref)
```


<br>

***

### Detection of Somatic Copy Number Alterations from Single-Cell Multiomics Data

The identification of somatic copy number alterations (CNAs) in cancer cells is crucial for understanding tumor evolution, including clonal dynamics causing relapse, and identifying potential therapeutic targets. While existing tools provide valuable insights into subclonal CNAs, they are typically limited to analyzing one type of omics data. In response to the growing use of cutting-edge technologies enabling simultaneous sequencing of multiple omics from individual cells, there emerges a need for new approaches that leverage multiomics data integration to improve the detection of CNAs. 

Addressing this need, we developed an R package, *muscadet*, that integrates multiple single-cell datasets across different omics modalities to enhance the accuracy and resolution of CNA detection within tumoral subclones. We demonstrated the potency of our approach through the analysis of single-cell Multiome data, integrating both single-cell RNA-seq and single-cell ATAC-seq datasets from a common pool of cells, across several multiple myeloma samples. *muscadet* outperformed existing copy number analysis tools for both scRNA-seq and scATAC-seq data, revealing accurate CNA profiles and subclones, validated by matched whole genome sequencing data. 

By providing a unified CNA analysis framework applicable to any combination of single-cell omics data, *muscadet* empowers researchers to unravel the clonal structure of tumor samples and uncover complex genomic alterations driving cancer progression.

***
