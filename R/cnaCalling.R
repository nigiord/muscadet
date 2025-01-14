

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

    procSnps <- utils::getFromNamespace("procSnps", "facets")

    # Step 1: Process SNPs and retain cluster & signal columns
    pmat_1 <- procSnps(
        rcmat,
        ndepth = ndepth,
        het.thresh = het.thresh,
        snp.nbhd = snp.nbhd,
        nX = nX,
        unmatched = FALSE,
        ndepthmax = ndepthmax
    )
    colnames(pmat_1)[7:8] <- c("cluster", "signal")  # Rename columns

    # Step 2: Process SNPs with only required columns to get het & keep columns
    pmat_2 <- procSnps(
        rcmat[, 1:6],
        ndepth = ndepth,
        het.thresh = het.thresh,
        snp.nbhd = snp.nbhd,
        nX = nX,
        unmatched = FALSE,
        ndepthmax = ndepthmax
    )

    # Combine both results
    pmat <- cbind(pmat_1, pmat_2[, c("het", "keep")])

    # Step 3: Compute logR and logOR with GC correction for logR
    counts2logROR <- utils::getFromNamespace("counts2logROR", "facets")
    dmat <- counts2logROR(pmat[pmat$rCountT > 0, ], gbuild, unmatched = FALSE)

    # Step 4: Exclude log ratios for allelic data, replace by mean of neighbors
    dmat <- dmat %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$chrom) %>% # Group by chromosome
        dplyr::arrange(.data$maploc, .by_group = TRUE) %>% # Ensure rows are ordered by maploc
        dplyr::mutate(
            cnlr = dplyr::if_else(
                .data$signal == "allelic",
                # Compute mean of neighbors
                (dplyr::lag(.data$cnlr, default = 0) + dplyr::lead(.data$cnlr, default = 0)) /
                    (ifelse(is.na(lag(.data$cnlr, default = NA)), 0, 1) + ifelse(is.na(lead(.data$cnlr, default = NA)), 0, 1)),
                .data$cnlr # Keep existing values for coverage rows
            )
        ) %>%
        dplyr::ungroup()
    dmat <- as.data.frame(dmat)

    # Step 5: Segment SNPs
    segsnps <- utils::getFromNamespace("segsnps", "facets")
    seg_results <- segsnps(dmat,
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


#' Get consensus segments across clusters
#'
#' This function processes segments data, in which clusters have different
#' segment breakpoints, to identify consensus segments across all clusters. It
#' groups breakpoints within a specified genomic distance (`dist.breakpoints`)
#' and calculates representative breakpoints for each group based on cluster
#' size and median coordinates.
#'
#' @param x A data frame containing cluster segments information with the
#'   following required columns:
#'   - `chrom`: Chromosome name (`factor` or `character`).
#'   - `start`: Start position of the segment (`numeric`).
#'   - `end`: End position of the segment (`numeric`).
#'   - `cluster`: Cluster identifier for each breakpoint (`numeric`).
#' @param ncells A named vector specifying the number of cells per cluster
#'   (`numeric` vector). The names should correspond to the cluster identifiers
#'   in the `cluster` column of `x`.
#' @param dist.breakpoints A numeric value specifying the minimum genomic
#'   distance between adjacent breakpoints to be grouped into the same consensus
#'   segment (`numeric` value). Default: `10e6`.
#'
#' @return A data frame containing consensus segments with the following columns:
#'   - `chrom`: Chromosome name.
#'   - `start`: Start position of the consensus segment.
#'   - `end`: End position of the consensus segment.
#'
#' @import dplyr
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges reduce
#' @importFrom S4Vectors subjectHits
#' @importFrom rlang .data
#'
#' @export
#'
#' @seealso [muscadet::cnaCalling()]
#'
#' @examples
#' # Example data frame
#' segs <- data.frame(
#'   chrom = c("chr1", "chr1", "chr2", "chr2"),
#'   start = c(1.2e6, 1.1e6, 3.1e6, 3.2e6),
#'   end = c(2.5e6, 2.6e6, 5.5e6, 5.7e6),
#'   cluster = c("1", "2", "1", "2")
#' )
#'
#' # Generate consensus segments
#' consensus_segs <- getSegConsensus(segs,
#'                                   ncells = c("1" = 50, "2" = 30),
#'                                   dist.breakpoints = 1e6)
#' print(consensus_segs)
#'
getSegConsensus <- function(x, ncells, dist.breakpoints = 10e6) {
    # Check if x is a data frame
    if (!is.data.frame(x)) {
        stop("'x' must be a data frame.")
    }

    # Check for required columns in x
    required_columns <- c("chrom", "start", "end", "cluster")
    missing_columns <- setdiff(required_columns, colnames(x))
    if (length(missing_columns) > 0) {
        stop(
            "The data frame 'x' is missing the following required columns: ",
            paste(missing_columns, collapse = ", ")
        )
    }

    # Check if dist.breakpoints is numeric and positive
    if (!is.numeric(dist.breakpoints) ||
        length(dist.breakpoints) != 1 || dist.breakpoints <= 0) {
        stop("'dist.breakpoints' must be a single positive numeric value.")
    }

    # Check for missing values in required columns
    if (anyNA(x[, required_columns])) {
        stop("The data frame 'x' contains missing values in required columns.")
    }

    # Get start and end breakpoints
    x$chrom <- factor(x$chrom, levels = unique(x$chrom))
    breakpoints_start <- dplyr::arrange(x, .data$chrom, .data$start)[, c("chrom", "start", "start", "cluster")]
    colnames(breakpoints_start)[c(2, 3)] <- c("start", "end")
    breakpoints_end <- dplyr::arrange(x, .data$chrom, .data$end)[, c("chrom", "end", "end", "cluster")]
    colnames(breakpoints_end)[c(2, 3)] <- c("start", "end")

    # Merge all breakpoints in one GRanges object
    breakpoints <- GenomicRanges::GRanges(arrange(
        rbind(breakpoints_start, breakpoints_end),
        .data$chrom,
        .data$start,
        .data$end
    ))

    # Find clusters of breakpoints (named "group" here to distinguish from the clusters of cells)
    breakpoints$group <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(
        breakpoints,
        GenomicRanges::reduce(breakpoints, min.gapwidth = dist.breakpoints)
    ))
    breakpoints <- as.data.frame(breakpoints)

    # Order clusters of cells from the smallest to the largest
    breakpoints$cluster <- factor(breakpoints$cluster,
                                  levels = names(ncells)[order(ncells)],
                                  ordered = T)

    # In case of unique breakpoint in a chromosome (unique segment < dist.breakpoints):
    # separation into 2 distinct new groups
    breakpoints <- breakpoints %>%
        dplyr::group_by(.data$seqnames, .data$cluster) %>%
        dplyr::mutate(ngroupperchr = length(unique(.data$group))) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$seqnames, .data$start) %>%
        dplyr::mutate(group = case_when(
            .data$ngroupperchr == 1 ~ max(breakpoints$group) + cur_group_id(),
            .default = .data$group
        )) %>%
        dplyr::ungroup()

    # If there are several breakpoints in one group of breakpoints per cluster:
    # compute median to get a unique breakpoint
    breakpoints <- breakpoints %>%
        dplyr::group_by(.data$group, .data$cluster) %>%
        dplyr::mutate(start = median(.data$start), end = median(.data$end)) %>%
        dplyr::ungroup() %>%
        unique()

    # In each group of breakpoints: select the breakpoint coordinates from the largest cluster of cells
    breakpoints <- breakpoints %>%
        dplyr::group_by(.data$group) %>%
        dplyr::mutate(start = .data$start[.data$cluster == max(.data$cluster)],
                      end = .data$end[.data$cluster == max(.data$cluster)]) %>%
        dplyr::select(-.data$cluster, -.data$ngroupperchr) %>%
        unique()

    # Transform breakpoints back into segments
    colnames(breakpoints)[1] <- "chrom"
    consensus_segs <- breakpoints %>%
        dplyr::group_by(.data$chrom) %>%
        dplyr::mutate(start = if (n() == 1) {
            .data$start # Single breakpoint: start remains as is
        } else {
            c(.data$start[1:(length(.data$start) - 1)], NA) # Multiple breakpoints: create segments
        }, end = if (n() == 1) {
            .data$end # Single breakpoint: end remains as is
        } else {
            c(.data$end[2:length(.data$end)], NA) # Multiple breakpoints: create segments
        }) %>%
        dplyr::select(-.data$group, -.data$width, -.data$strand) %>%
        dplyr::filter(!is.na(.data$start) & !is.na(.data$end)) # Remove rows with NA start or end

    return(as.data.frame(consensus_segs))
}

