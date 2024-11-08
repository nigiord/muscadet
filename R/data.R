# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Cell barcodes ----------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Cell barcodes
#'
#' @description
#' Character vectors of unique cell barcodes.
#'
#' @name barcodes
#' @rdname barcodes
#'
#' @format
#' A character vector.
#'
"barcodes_atac_tumor"

#' @name barcodes
#' @rdname barcodes
#' @format NULL
#'
"barcodes_atac_ref"

#' @name barcodes
#' @rdname barcodes
#' @format NULL
#'
"barcodes_rna_tumor"

#' @name barcodes
#' @rdname barcodes
#' @format NULL
#'
"barcodes_rna_ref"

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Feature coordinates ----------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Feature coordinates
#'
#' @description
#' Data frames of features (peaks, genes...) coordinates on genome.
#'
#' @name features
#' @rdname features
#'
#' @format
#' A data frame with the following columns:
#' \describe{
#'   \item{`CHROM`}{Chromosome in chracter format, e.g. "15", "X" (`character`).}
#'   \item{`start`}{Start position of the feature (`integer`).}
#'   \item{`end`}{End position of the feature (`character`).}
#'   \item{`id`}{Unique identifier, e.g. gene name "CDH1" or peak identifier
#'   CHROM_start_end "1_1600338_1600838" (`character`).}
#' }
#'
"genes"

#' @name features
#' @rdname features
#' @format NULL
#'
"peaks"


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Matrices of counts -----------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Matrices of raw counts
#'
#' @description
#' Matrices of raw read counts *features x cells* in `dgCMatrix` format.
#'
#' @name mat_counts
#' @rdname mat_counts
#'
#' @format
#' A `dgCMatrix` matrix of numeric values with the following dimensions:
#' \describe{
#'   \item{`rows`}{Features (peaks, genes).}
#'   \item{`columns`}{Cell barcodes.}
#' }
#'
"mat_counts_atac_tumor"

#' @name mat_counts
#' @rdname mat_counts
#' @format NULL
#'
"mat_counts_atac_ref"

#' @name mat_counts
#' @rdname mat_counts
#' @format NULL
#'
"mat_counts_rna_tumor"

#' @name mat_counts
#' @rdname mat_counts
#' @format NULL
#'
"mat_counts_rna_ref"

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Allele counts ----------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Allele counts at variation positions
#'
#' @description
#' Data frames of allele counts at single nucleotide polymorphisms (SNPs)
#' positions per cell.
#'
#' @name allele_counts
#' @rdname allele_counts
#'
#' @format
#' A data frame with columns based on the [variant calling format (VCF)
#' format](https://en.wikipedia.org/wiki/Variant_Call_Format) columns.
#' It contains the following columns:
#' \describe{
#'   \item{`cell`}{Barcodes of cells (`character`).}
#'   \item{`id`}{SNP unique identifier defined as
#'   CHROM_POS_REF_ALT, e.g. "1_920949_C_G" (`character`).}
#'   \item{`CHROM`}{Chromosome in integer format, e.g. 15 (X and Y chromosomes
#'   are not included) (`integer`).}
#'   \item{`POS`}{Position of the SNP (1-base positions) (`integer`).}
#'   \item{`REF`}{Reference allele base, "A" "C" "G" or "T" (`character`).}
#'   \item{`ALT`}{Alternative allele base, "A" "C" "G" or "T" (`character`).}
#'   \item{`RD`}{Reference allele depth/count (`integer`).}
#'   \item{`AD`}{Alternative allele depth/count (`integer`).}
#'   \item{`DP`}{Total depth/count (`integer`).}
#'   \item{`GT`}{Genotype, "0|1" or 1|0" (can be "0/1" if unphased) (`character`).}
#' }
#'
"allele_counts_atac_tumor"

#' @name allele_counts
#' @rdname allele_counts
#' @format NULL
#'
"allele_counts_atac_ref"

#' @name allele_counts
#' @rdname allele_counts
#' @format NULL
#'
"allele_counts_rna_tumor"

#' @name allele_counts
#' @rdname allele_counts
#' @format NULL
#'
"allele_counts_rna_ref"


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Log ratio from bulk data -----------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Log R ratio from bulk sequencing data
#'
#' @description Data frame containing log R ratio values per genomic segments
#' from bulk sequencing data.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{`CHROM`}{Chromosome in integer format, e.g. 15, 23 (for X chromosome) (`character`).}
#'   \item{`start`}{Start position of the segment (`integer`).}
#'   \item{`end`}{End position of the segment (`character`).}
#'   \item{`lrr`}{Log R ratio of the segment ("cnlr.median" column from
#'   [facets::fitcncf()] `$cncf` data frame) (`numeric`).}
#' }
#'
#' @note Data obtained from whole genome sequencing (WGS) after using
#' [facets::fitcncf()] from [facets] "Cellular Fraction and Copy Numbers from
#' Tumor Sequencing" version `0.6.2`: `$cncf` data frame columns `chrom`, `start`,
#' `end`, and `cnlr.median`.
#'
#' @references
#' \describe{
#'   \item{[facets] package}{Shen R, Seshan VE. FACETS: allele-specific copy number and
#'   clonal heterogeneity analysis tool for high-throughput DNA sequencing.
#'   Nucleic Acids Res. 2016 Sep 19;44(16):e131.
#'   doi: [10.1093/nar/gkw520](https://www.doi.org/10.1093/nar/gkw520).
#'   PMID: 27270079; PMCID: PMC5027494.}
#' }
#'
"bulk_lrr"


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Genome chromosome sizes -----------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Genome chromosome sizes
#'
#' @description \code{\link{GRanges}} objects containing chromosomes sizes for
#'   hg38, hg19 and mm10 genome assemblies.
#'
#' @name genome_chrom
#' @rdname genome_chrom
#'
#' @format A \code{\link{GRanges}} object containing:
#' \describe{
#'   \item{`seqnames`}{Chromosome name: 1 to 22, X and Y for human ; and 1 to
#'   19, X and Y for mouse (`Rle`).}
#'   \item{`ranges`}{Ranges of chromosomes (`IRanges`).}
#'   \item{`strand`}{Strand information * (`Rle`).}
#' }
#'
#' @note Data obtained from assemblies provided by the [BSgenome] package.
#' - `BSgenome.Hsapiens.UCSC.hg38` version 1.4.5 - `GRCh38.p14`
#' - `BSgenome.Hsapiens.UCSC.hg19` version 1.4.3 - `GRCh37.p13`
#' - `BSgenome.Mmusculus.UCSC.mm10` version 1.4.3 - `GRCm38.p6`
#'
"hg38_chrom"

#' @name genome_chrom
#' @rdname genome_chrom
#' @format NULL
#'
"hg19_chrom"

#' @name genome_chrom
#' @rdname genome_chrom
#' @format NULL
#'
"mm10_chrom"
