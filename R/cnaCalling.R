

#' Process read count matrix and segmentation
#'
#' Function adapted from [facets::preProcSample()] to process read count matrix
#' and generates a segmentation tree. The modifications from the original
#' function includes:
#' - Incorporation of `cluster` and `signal` columns into the final result.
#' - Change in the log ratio (`cnlr`) for allelic data, it is computed as the
#' log ratio (`cnlr`) mean between the previous and next coverage positions.
#'
#' @param rcmat A data frame with 8 required columns: `Chrom`, `Pos`, `NOR.DP`,
#'   `NOR.RD`, `TUM.DP`, `TUM.RD`, `cluster`, and `signal` (`data.frame`).
#' @param het.thresh VAF (Variant Allele Frequency) threshold to call SNP
#'   heterozygous (`numeric`). Default: 0.25.
#' @param snp.nbhd Window size for selecting loci to reduce serial correlation
#'   (`numeric`). Default: 250.
#' @param cval Critical value for segmentation (`numeric`). Default: 25.
#' @param gbuild Genome build used for alignment. One of \code{"hg19"},
#'   \code{"hg38"}, or \code{"mm10"} (`character`). Default: \code{"hg19"}.
#' @param hetscale Logical value indicating whether log odds ratio (logOR)
#'   should be scaled to give more weight in the test statistics for
#'   segmentation and clustering (`logical`). Usually only 10% of snps are hets
#'   and hetscale gives the logOR contribution to T-square as 0.25/proportion of
#'   hets. Default: \code{TRUE}.
#' @param ndepth Minimum depth in normal reference to keep (`numeric`). Default:
#'   5.
#' @param ndepthmax Maximum normal coverage threshold for filtering loci
#'   (`numeric`). Default: 1000.
#'
#' @return A list containing:
#' \describe{
#'   \item{pmat}{Read counts and other elements of all the loci.}
#'   \item{gbuild}{Genome build used for the analysis.}
#'   \item{nX}{Chromosome number for X (e.g., 23 for human, 20 for mouse).}
#'   \item{clusters}{Unique clusters from the processed data.}
#'   \item{seg.tree}{Segmentation tree for each chromosome.}
#'   \item{jointseg}{Segmented SNP data.}
#'   \item{hscl}{Scaling factor for logOR data.}
#' }
#' @details
#' The function processes SNP data to generate a segmentation tree. It uses
#' \code{\link[facets]{procSnps}} to compute initial values, adjusts the log ratio
#' for allelic signals, and computes segmentation using \code{\link[facets]{segsnps}}.
#'
#' SNPs in a genome are not evenly spaced, and loci are sampled within the specified
#' window (\code{snp.nbhd}) to reduce serial correlation.
#'
#' @seealso \code{\link[facets]{preProcSample}}, \code{\link[facets]{procSnps}},
#'   \code{\link[facets]{segsnps}}
#'
#' @source This function is derived from the [facets::preProcSample()] function,
#' with the modifications to fit its use in muscadet.
#'
#' Seshan VE, Shen R (2021). _facets: Cellular Fraction and Copy Numbers from Tumor Sequencing_.
#' R package version 0.6.2, [https://github.com/mskcc/facets](https://github.com/mskcc/facets).
#'
#' @references Shen R, Seshan VE. FACETS: allele-specific copy number and clonal
#' heterogeneity analysis tool for high-throughput DNA sequencing. Nucleic Acids
#' Res. 2016 Sep 19;44(16):e131. doi: [10.1093/nar/gkw520](http://doi.org/10.1093/nar/gkw520).
#'
#' @import dplyr
#' @import facets
#'
#' @export
#'
#' @examples
#' # Load muscadet object
#' data(muscadet_obj)
#' counts <- muscadet_obj$cnacalling$combined.counts
#' counts <- counts[complete.cases(counts),]
#' counts_clus <- counts[which(counts$cluster == 1),]
#' result <- preProcSample2(counts_clus)
#'
preProcSample2 <- function(
        rcmat,
        het.thresh = 0.25,
        snp.nbhd = 250,
        cval = 25,
        gbuild = "hg38",
        hetscale = TRUE,
        ndepth = 5,
        ndepthmax = 1000
) {

    # rcmat correct format
    rcmat <- as.data.frame(rcmat)

    gbuild <- match.arg(gbuild, c("hg19", "hg38", "mm10"))

    # Determine chromosome number for X based on genome build
    nX <- switch(
        gbuild,
        "hg19" = 23,
        "hg38" = 23,
        "mm10" = 20,
        stop("Unsupported genome: ", gbuild)
    )

    # Step 1: Process SNPs and retain cluster & signal columns
    pmat_1 <- facets:::procSnps(rcmat,
                                ndepth = ndepth,
                                het.thresh = het.thresh,
                                snp.nbhd = snp.nbhd,
                                nX = nX,
                                unmatched = FALSE,
                                ndepthmax = ndepthmax)
    colnames(pmat_1)[7:8] <- c("cluster", "signal")  # Rename columns

    # Step 2: Process SNPs with only required columns to get het & keep columns
    pmat_2 <- facets:::procSnps(rcmat[, 1:6],
                                ndepth = ndepth,
                                het.thresh = het.thresh,
                                snp.nbhd = snp.nbhd,
                                nX = nX,
                                unmatched = FALSE,
                                ndepthmax = ndepthmax)

    # Combine both results
    pmat <- cbind(pmat_1, pmat_2[, c("het", "keep")])

    # Step 3: Compute logR and logOR with GC correction for logR
    dmat <- facets:::counts2logROR(pmat[pmat$rCountT > 0, ], gbuild, unmatched)

    # Step 4: Exclude log ratios for allelic data, replace by mean of neighbors
    dmat <- dmat %>%
        dplyr::ungroup() %>%
        dplyr::group_by(chrom) %>% # Group by chromosome
        dplyr::arrange(maploc, .by_group = TRUE) %>% # Ensure rows are ordered by maploc
        dplyr::mutate(
            cnlr = dplyr::if_else(
                signal == "allelic",
                # Compute mean of neighbors
                (dplyr::lag(cnlr, default = 0) + dplyr::lead(cnlr, default = 0)) /
                    (ifelse(is.na(lag(cnlr, default = NA)), 0, 1) + ifelse(is.na(lead(cnlr, default = NA)), 0, 1)),
                cnlr # Keep existing values for coverage rows
            )
        ) %>%
        dplyr::ungroup()
    dmat <- as.data.frame(dmat)

    # Step 5: Segment SNPs
    seg_results <- facets:::segsnps(dmat,
                                    cval = cval,
                                    hetscale = hetscale,
                                    delta = 0)

    # Make sure the joint segmentation data frame is ordered
    seg_results$jointseg <- seg_results$jointseg[order(seg_results$jointseg[, "maploc"]), ]
    seg_results$jointseg <- seg_results$jointseg[order(seg_results$jointseg[, "chrom"]), ]

    # Prepare output
    result <- list(
        pmat = pmat,
        gbuild = gbuild,
        nX = nX,
        clusters = unique(pmat$cluster)
    )
    if(!is.null(pmat$cluster)) result$clusters <- unique(pmat$cluster)

    # Combine output with segmentation results
    return(c(result, seg_results))
}

