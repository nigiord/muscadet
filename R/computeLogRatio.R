#' Compute log R ratios
#'
#' Computes log R ratios from raw count matrices within \code{\link{muscadet}}
#' objects. The log R ratios values are computed based on read counts from a
#' sample \code{\link{muscadet}} object *versus* read counts a the `reference`
#' \code{\link{muscadet}} object. In the output sample \code{\link{muscadet}}
#' object, the newly computed matrix of log R ratios is added for the
#' selected `omic`.
#'
#' @param x A \code{\link{muscadet}} object containing sample data (`muscadet`).
#' @param reference Another \code{\link{muscadet}} object containing reference
#'   data (`muscadet`).
#' @param omic Name of the omic to apply this function (`character` string).
#' @param method Method to apply to the selected omic (`character` string).
#'   Supported methods are "ATAC" and "RNA":
#'   * "ATAC" method calls for [computeLogRatioATAC()] function
#'   * "RNA" method calls for [computeLogRatioRNA()] function
#'   If `NULL` and if the omic type is either "ATAC" or "RNA", the corresponding
#'   method will be applied.
#' @param new.label.features New label for features (`character` string). If
#'   `NULL`, the label remains unchanged when using the "RNA" method and becomes
#'   "windows of peaks" when using the "ATAC" method.
#' @param remove.raw `TRUE` or `FALSE` (`logical`). Whether to remove raw count
#'   matrices. `TRUE` by default to reduce object size. Setting it to `FALSE`
#'   will keep raw count matrices within the object.
#' @param quiet `TRUE` or `FALSE` (`logical`). Whether to turn off messages. By
#'   default: `FALSE`.
#'
#' @inheritDotParams computeLogRatioATAC windowSize slidingSize minReads
#'   minPeaks thresh_capping
#' @inheritDotParams computeLogRatioRNA genesPerWindow refReads refMeanReads
#'
#' @return
#' A \code{\link{muscadet}} object corresponding to the sample
#' \code{\link{muscadet}} object (`x`) containing the computed log R ratio
#' matrix in the `coverage` slot of the selected `omic`.
#'
#' @details
#' Log R ratios computation steps are described in [computeLogRatioATAC()]
#' and [computeLogRatioRNA()] functions.
#'
#' @family computeLogRatio
#'
#' @export
#'
#' @examples
#' # compute log R ratios for ATAC
#' muscadet <- computeLogRatio(
#'   x = muscadet,
#'   reference = muscadet_ref,
#'   omic = "ATAC",
#'   method = "ATAC",
#'   minReads = 1, # low value only for the example subsampled datasets
#'   minPeaks = 10 # low value only for the example subsampled datasets
#' )
#'
#' # compute log R ratios for RNA
#' muscadet <- computeLogRatio(
#'   x = muscadet,
#'   reference = muscadet_ref,
#'   omic = "RNA",
#'   method = "RNA",
#'   refReads = 20
#' )
#'
computeLogRatio <- function(x,
                            reference,
                            omic,
                            method = NULL,
                            new.label.features = c(
                              ATAC = "windows of peaks",
                              RNA = "genes"
                            ),
                            remove.raw = TRUE,
                            quiet = FALSE,
                            ...) {

    # Check omic is in the omics of 'x' and 'reference'
    stopifnot(
        "'omic' argument must corresponds to an omic name in both muscadet objects 'x' and 'reference'." =
            omic %in% names(x@omics) & omic %in% names(reference@omics)
    )

  # Check method
  if (is.null(method)) method <- x@omics[[omic]]@type
  stopifnot(
    "'method' must be either 'ATAC' or 'RNA'." =
      method %in% c("ATAC", "RNA")
  )

  # Check raw count matrix in muscadet object
  stopifnot(
      "Raw count matrix not found in the muscadet 'x' object." =
          !is.null(x@omics[[omic]]@coverage[["mat.counts"]])
  )


  # ATAC method ----------------------------------------------------------------

  if (method == "ATAC") {

      # Check dot arguments
      dots <- list(...)
      stopifnot("Unknown arguments used." =
                    all(names(dots) %in% c(formalArgs(computeLogRatioATAC))))

      # Check label features
      if (is.null(new.label.features)) {
          new.label.features <- c(ATAC = "windows of peaks")
      }

      if (quiet == FALSE) {
          message("-- computeLogRatio: Method 'ATAC' using computeLogRatioATAC().")
      }

      obj <- computeLogRatioATAC(
          matTumor = x@omics[[omic]]@coverage[["mat.counts"]],
          matRef = reference@omics[[omic]]@coverage[["mat.counts"]],
          peaksCoord = x@omics[[omic]]@coverage[["coord.features"]],
          quiet = quiet,
          ...
      )
  }

  # RNA method ---------------------------------------------------------------

  if (method == "RNA") {

    # Check dot arguments
    dots <- list(...)
    stopifnot(
      "Unknown arguments used." =
        all(names(dots) %in% c(formalArgs(
          computeLogRatioRNA
        )))
    )

    # Check label features
    if (is.null(new.label.features)) {
      new.label.features <- c(RNA = x@omics[[omic]]@coverage[["label.features"]])
    }

    if (quiet == FALSE) {
        message("-- computeLogRatio: Method 'RNA' using computeLogRatioRNA().")
    }

    obj <- computeLogRatioRNA(
      matTumor = x@omics[[omic]]@coverage[["mat.counts"]],
      matRef = reference@omics[[omic]]@coverage[["mat.counts"]],
      genesCoord = x@omics[[omic]]@coverage[["coord.features"]],
      quiet = quiet,
      ...
    )
  }

  if (exists("all_steps")) {
      if (all_steps == TRUE) {
          obj <- list(
              matTumor = obj$step08$matTumor,
              matRef = obj$step08$matRef,
              params = obj$params,
              coord = obj$coord
          )
      }
  }

  x@omics[[omic]]@coverage[["log.ratio"]] <- obj$matTumor
  x@omics[[omic]]@coverage[["coord.features"]] <- obj$coord
  x@omics[[omic]]@coverage[["label.features"]] <- new.label.features[omic]

  if (remove.raw == TRUE) {
    x@omics[[omic]]@coverage[["mat.counts"]] <- NULL
  }

  if (quiet == FALSE) {
      message("Done.")
  }

  return(x)
}


#' Compute log R ratios for scATAC-seq data
#'
#' Compute log R ratios from raw count matrices with method specifically adapted
#' for scATAC-seq data. The counts per peaks are grouped into counts per windows
#' of peaks, thereby features become windows of peaks instead of peaks.
#'
#' @param matTumor Raw count matrix *features x cells* for tumor/sample cells
#'   (`matrix` or `dgCMatrix`).
#' @param matRef Raw count matrix *features x cells* for reference cells
#'   (`matrix` or `dgCMatrix`).
#' @param peaksCoord Data frame of peak coordinates with columns `CHROM`, `start`,
#'   `end`, `id` (`data.frame`).
#' @param genome Reference genome name among: "hg38", "hg19" and "mm10"
#'   (`character` string). By default: "hg38".
#' @param windowSize Size of windows in base pairs (`integer` value). By default: `10e6` (10
#'   Mbp).
#' @param slidingSize Distance between start positions of sliding windows in base pairs
#'   (`integer` value). If set to the same value as `windowSize`, the windows
#'   don't overlap. By default: `2e6` (2 Mbp).
#' @param minReads Minimum read average per window (`integer` value). By
#'   default: `5`.
#' @param minPeaks Minimum number of peaks per window (`integer` value). By
#'   default: `100`.
#' @param thresh_capping Threshold to cap the range of log R ratio values
#'   (`numeric` value). By default: `3`.
#' @param all_steps `TRUE` or `FALSE` (`logical`). Whether to keep intermediate
#'   result from every step in the final object. By default: `FALSE`.
#' @param quiet `TRUE` or `FALSE` (`logical`). Whether to turn off messages. By
#'   default: `FALSE`.
#'
#'
#' @return
#' If `all_steps` is set to `FALSE`, a list containing:
#' \describe{
#'  \item{`matTumor`}{Matrix of log R ratio values *features x cells* for the
#'  tumor/sample cells (`matrix`).}
#'  \item{`matRef`}{Matrix of log R ratio values *features x cells* for the
#'  reference cells (`matrix`).}
#'  \item{`params`}{List of parameters set for the `windowSize`, `slidingSize`,
#'  `minReads`, `minPeaks` and `thresh_capping` arguments (`list`).}
#'  \item{`coord`}{Data frame of coordinates for windows of peaks and associated
#'  data along the different steps (`data.frame`).
#'  Columns :
#'  - `CHR`, `start`, `end`, `width`, `id`: coordinates and unique identifier of
#'  windows (depends on `windowSize` and `slidingSize` arguments).
#'  - `nPeaks`: number of peaks per window.
#'  - `sumReads.tum/ref`: sum of read counts for all cells in tumor or reference
#'  cells.
#'  - `meanReads.tum/ref`: mean of read counts per cells for tumor or reference
#'  cells.
#'  - `sdReads.tum/ref`: standard deviation of read counts per cells for tumor
#'  or reference cells.
#'  - `keep`: logical, `TRUE` for windows to keep after filtering based on
#'  coverage (depends on `minPeaks` and `meanReads` arguments).
#'  - `meanReads/sdReads.norm.tum/ref`: mean/sd of normalized counts per million
#'  for tumor/reference cells.
#'  - `meanLRR/sdReads.raw.tum/ref`: mean/sd of raw log R ratio (LRR) for
#'  tumor/reference cells.
#'  - `meanLRR/sdLRR.cap.tum/ref`: mean/sd of capped log R ratio (LRR) for
#'  tumor/reference cells (depends on `thresh_capping` argument).
#'  - `meanLRR/sdLRR.cent.tum/ref`: mean/sd of centered log R ratio (LRR) for
#'  tumor/reference cells.
#'  - `meanLRR/sdLRR.corr.tum/ref`: mean/sd of final log R ratio (LRR) corrected
#'  by reference variability for tumor/reference cells.}
#' }
#' If `all_steps` is set to `TRUE`, the previously described list can be found for
#' each step in a different list element.
#'
#' @family computeLogRatio
#'
#' @import Matrix
#' @import dplyr
#' @importFrom BiocGenerics width
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges slidingWindows
#' @importFrom GenomicRanges findOverlaps
#'
#' @export
#'
#' @details
#' The raw count matrix is transformed into log R ratios through the following
#' steps:
#'  - Group peaks in windows (`windowSize` and `slidingSize` arguments)
#'  - Filtering on coverage (`minPeaks` and `meanReads` arguments)
#'  - Normalization for sequencing depth
#'  - Log transformation and normalization by reference data: log R ratio
#'  - Capping the range of values (`thresh_capping` argument)
#'  - Centering of cells
#'  - Correcting by reference variability
#'
#' @examples
#' obj_atac <- computeLogRatioATAC(
#'   matTumor = matCounts(muscadet)$ATAC,
#'   matRef = matCounts(muscadet_ref)$ATAC,
#'   peaksCoord = coordFeatures(muscadet)$ATAC,
#'   genome = slot(muscadet, "genome"),
#'   minReads = 1, # low value only for the example subsampled datasets
#'   minPeaks = 10 # low value only for the example subsampled datasets
#' )
#' table(obj_atac$coord$keep)
#'
#' # With results form every step when `all_steps = TRUE`
#' obj_atac_all <- computeLogRatioATAC(
#'   matTumor = matCounts(muscadet)$ATAC,
#'   matRef = matCounts(muscadet_ref)$ATAC,
#'   peaksCoord = coordFeatures(muscadet)$ATAC,
#'   genome = slot(muscadet, "genome"),
#'   minReads = 1, # low value only for the example subsampled datasets
#'   minPeaks = 10, # low value only for the example subsampled datasets
#'   all_steps = TRUE
#' )
#' names(obj_atac_all)
#' table(obj_atac_all$coord$keep)
#' nrow(obj_atac_all$step08$matTumor)
#'
computeLogRatioATAC <- function(matTumor,
                                matRef,
                                peaksCoord,
                                genome = "hg38",
                                windowSize = 10e6,
                                slidingSize = 2e6,
                                minReads = 5,
                                minPeaks = 100,
                                thresh_capping = 3,
                                all_steps = FALSE,
                                quiet = FALSE) {
  params <- list(
    windowSize = windowSize,
    slidingSize = slidingSize,
    minReads = minReads,
    minPeaks = minPeaks,
    thresh_capping = thresh_capping
  )

  # Sort cells (column names of matrices)
  matTumor <- matTumor[, sort(colnames(matTumor))]
  matRef <- matRef[, sort(colnames(matRef))]

  ## Step 01: Group peaks in windows -------------------------------------------
  if (quiet == FALSE) {
    message(
      "Step 01 - Group peaks in windows: window size set at ",
      windowSize / 1000000,
      " Mb, sliding by ",
      slidingSize / 1000000,
      " Mb"
    )
  }

  # Select genome
  if (genome == "hg38") {
    genome_chrom <- muscadet:::hg38_chrom
  }
  if (genome == "hg19") {
    genome_chrom <- muscadet:::hg19_chrom
  }
  if (genome == "mm10") {
    genome_chrom <- muscadet:::mm10_chrom
  }

  # Create windows
  windowsGR <- GenomicRanges::slidingWindows(
    x = genome_chrom,
    width = windowSize,
    step = slidingSize
  ) %>%
    unlist() %>%
    .[which(BiocGenerics::width(.) >= windowSize / 10), ] # keep windows at the end of chr that are at least 1/10 of window width

  # Match peaks and windows
  peaksCoord <- GenomicRanges::GRanges(peaksCoord)
  ovlaps <- as.data.frame(GenomicRanges::findOverlaps(windowsGR, peaksCoord, ignore.strand = TRUE))

  # Get the sum of read counts per window
  window_index <- unique(ovlaps[, 1])
  matTumor <- t(sapply(window_index, function(i) {
    peak_index <- ovlaps[ovlaps[, 1] == i, 2]
    mat <- matTumor[peak_index, ]
    if (is.vector(mat)) {
      mat <- t(as.matrix(mat))
    } # case with only 1 peak in windows (vector instead of matrix)
    Matrix::colSums(mat)
  }))
  rownames(matTumor) <- window_index

  matRef <- t(sapply(window_index, function(i) {
    peak_index <- ovlaps[ovlaps[, 1] == i, 2]
    mat <- matRef[peak_index, ]
    if (is.vector(mat)) {
      mat <- t(as.matrix(mat))
    } # case with only 1 peak in windows (vector instead of matrix)
    Matrix::colSums(mat)
  }))
  rownames(matRef) <- window_index

  obj <- list()

  if (all_steps == T) {
    obj$step01 <- list(
      matTumor = matTumor,
      matRef = matRef,
      name = "Counts per windows"
    )
  }

  ## Step 02: Filtering on coverage --------------------------------------------
  if (quiet == FALSE) {
    message(
      "Step 02 - Filtering windows: Minimum of ",
      minPeaks,
      " peaks per window with a minimum average of ",
      minReads,
      " read(s)"
    )
  }

  # Initialize window coordinates data frame
  windowsDF <- as.data.frame(windowsGR)
  colnames(windowsDF)[1] <- "CHR"
  windowsDF[, "strand"] <- NULL

  # Add data to window coordinates data frame
  win <- 1:nrow(windowsDF) %in% window_index
  windowsDF$id <- 1:nrow(windowsDF)
  windowsDF$nPeaks[win] <- table(ovlaps[1])
  windowsDF$nPeaks[!win] <- 0
  windowsDF$sumReads.tum[win] <- Matrix::rowSums(matTumor)
  windowsDF$meanReads.tum[win] <- Matrix::rowMeans(matTumor)
  windowsDF$sdReads.tum[win] <- matrixStats::rowSds(matTumor)
  windowsDF$sumReads.ref[win] <- Matrix::rowSums(matRef)
  windowsDF$meanReads.ref[win] <- Matrix::rowMeans(matRef)
  windowsDF$sdReads.ref[win] <- matrixStats::rowSds(matRef)

  # Add filtering information 'keep' in window coordinates data frame
  windowsDF <- dplyr::mutate(windowsDF, keep = meanReads.ref >= minReads &
    windowsDF$nPeaks >= minPeaks)

  stopifnot("No windows passing filters (`minReads` and `minPeaks` arguments)" = any(windowsDF$keep))

  # Filter matrices
  matTumor <- matTumor[windowsDF$keep[window_index], ]
  matRef <- matRef[windowsDF$keep[window_index], ]

  if (all_steps == T) {
    obj$step02 <- list(
      matTumor = matTumor,
      matRef = matRef,
      name = "Counts per windows (filtered)"
    )
  }

  ## Step 03: Normalization for sequencing depth -------------------------------
  if (quiet == FALSE) {
    message("Step 03 - Normalization for sequencing depth: Normalized counts per million")
  }

  # Normalization and transformation in counts per million
  matTumor <- apply(matTumor, 2, function(x) {
    (x / sum(x, na.rm = T)) * 1e6
  })
  matRef <- apply(matRef, 2, function(x) {
    (x / sum(x, na.rm = T)) * 1e6
  })

  if (all_steps == T) {
    obj$step03 <- list(
      matTumor = matTumor,
      matRef = matRef,
      name = "Normalized counts per million"
    )
  }

  windowsDF$meanReads.norm.tum[windowsDF$keep] <- Matrix::rowMeans(matTumor)
  windowsDF$sdReads.norm.tum[windowsDF$keep] <- matrixStats::rowSds(matTumor)
  windowsDF$meanReads.norm.ref[windowsDF$keep] <- Matrix::rowMeans(matRef)
  windowsDF$sdReads.norm.ref[windowsDF$keep] <- matrixStats::rowSds(matRef)

  ## Step 04: Log transformation and normalization by reference data: log R ratio ----------
  if (quiet == FALSE) {
    message("Step 04 - Log transformation and normalization by reference data: log R ratio")
  }

  # Log transformation
  matTumor <- apply(matTumor, 2, function(x) {
    log2(1 + x)
  })
  matRef <- apply(matRef, 2, function(x) {
    log2(1 + x)
  })

  # Ratio versus the mean of reference cells
  matTumor <- apply(matTumor, 2, "-", Matrix::rowMeans(matRef))
  matRef <- apply(matRef, 2, "-", Matrix::rowMeans(matRef))

  if (all_steps == T) {
    obj$step04 <- list(
      matTumor = matTumor,
      matRef = matRef,
      name = "Log R ratio"
    )
  }

  windowsDF$meanLRR.raw.tum[windowsDF$keep] <- Matrix::rowMeans(matTumor)
  windowsDF$sdLRR.raw.tum[windowsDF$keep] <- matrixStats::rowSds(matTumor)
  windowsDF$meanLRR.raw.ref[windowsDF$keep] <- Matrix::rowMeans(matRef)
  windowsDF$sdLRR.raw.ref[windowsDF$keep] <- matrixStats::rowSds(matRef)


  ## Step 05: Capping the range of values --------------------------------------
  if (quiet == FALSE) {
    message("Step 05 - Capping the range of values: threshold = ", thresh_capping)
  }

  # Max and min thresholds of log R ratios
  matTumor[matTumor > thresh_capping] <- thresh_capping
  matTumor[matTumor < (-1 * thresh_capping)] <- -1 * thresh_capping
  matRef[matRef > thresh_capping] <- thresh_capping
  matRef[matRef < (-1 * thresh_capping)] <- -1 * thresh_capping

  if (all_steps == T) {
    obj$step05 <- list(
      matTumor = matTumor,
      matRef = matRef,
      name = "Log R ratio capped"
    )
  }

  windowsDF$meanLRR.cap.tum[windowsDF$keep] <- Matrix::rowMeans(matTumor)
  windowsDF$sdLRR.cap.tum[windowsDF$keep] <- matrixStats::rowSds(matTumor)
  windowsDF$meanLRR.cap.ref[windowsDF$keep] <- Matrix::rowMeans(matRef)
  windowsDF$sdLRR.cap.ref[windowsDF$keep] <- matrixStats::rowSds(matRef)

  ## Step 06 - [No step 06 for scATAC-seq] -------------------------------------
  if (quiet == FALSE) {
    message("Step 06 - [No step 06 for scATAC-seq]")
  }

  ## Step 07: Centering of cells -----------------------------------------------
  if (quiet == FALSE) {
    message("Step 07 - Centering of cells")
  }

  # Center cells
  meansTumor <- apply(matTumor, 2, function(x) {
    mean(x, na.rm = TRUE)
  })
  matTumor <- t(apply(matTumor, 1, "-", meansTumor))

  meansRef <- apply(matRef, 2, function(x) {
    mean(x, na.rm = TRUE)
  })
  matRef <- t(apply(matRef, 1, "-", meansRef))

  if (all_steps == T) {
    obj$step07 <- list(
      matTumor = matTumor,
      matRef = matRef,
      name = "Log R ratio centered"
    )
  }

  windowsDF$meanLRR.cent.tum[windowsDF$keep] <- Matrix::rowMeans(matTumor)
  windowsDF$sdLRR.cent.tum[windowsDF$keep] <- matrixStats::rowSds(matTumor)
  windowsDF$meanLRR.cent.ref[windowsDF$keep] <- Matrix::rowMeans(matRef)
  windowsDF$sdLRR.cent.ref[windowsDF$keep] <- matrixStats::rowSds(matRef)


  ## Step 08: Correcting by reference variability ------------------------------
  if (quiet == FALSE) {
    message("Step 08 - Correcting by reference variability")
  }

  sd.ref <- apply(matRef, 1, function(x) sd(x, na.rm = TRUE))
  mean.ref <- apply(matRef, 1, function(x) mean(x, na.rm = TRUE))

  matTumor <- apply(matTumor, 2, function(x) {
    (x - mean.ref) / sd.ref
  })
  matRef <- apply(matRef, 2, function(x) {
    (x - mean.ref) / sd.ref
  })

  if (all_steps == T) {
    obj$step08 <- list(
      matTumor = matTumor,
      matRef = matRef,
      name = "Log R ratio corrected by ref variablity"
    )
  }

  windowsDF$meanLRR.corr.tum[windowsDF$keep] <- Matrix::rowMeans(matTumor)
  windowsDF$sdLRR.corr.tum[windowsDF$keep] <- matrixStats::rowSds(matTumor)
  windowsDF$meanLRR.corr.ref[windowsDF$keep] <- Matrix::rowMeans(matRef)
  windowsDF$sdLRR.corr.ref[windowsDF$keep] <- matrixStats::rowSds(matRef)

  if (all_steps == T) {
    obj$params <- params
    obj$coord <- windowsDF
    return(obj)
  } else {
    obj_min <- list(
      matTumor = matTumor,
      matRef = matRef,
      params = params,
      coord = windowsDF
    )
    return(obj_min)
  }
}



#' Compute log R ratios for scRNA-seq data
#'
#' Compute log R ratios from raw count matrices with method specifically adapted
#' for scRNA-seq data.
#'
#' @param genesCoord Data frame of gene coordinates with columns `CHROM`,
#'   `start`, `end`, `id` (`data.frame`).
#' @param genome Reference genome name among: "hg38", "hg19" and "mm10"
#'   (`character`). By default: "hg38".
#' @param genesPerWindow Number of genes per moving window (`integer` value). By
#'   default: `101`.
#' @param refReads Minimum of reads in reference cells (`integer` value). By
#'   default: `100`.
#' @param refMeanReads Minimum of average reads per reference cell (`integer`
#'   value). By default: `0.01`.
#' @inheritParams computeLogRatioATAC
#'
#' @return
#' If `all_steps` is set to `FALSE`, a list containing:
#' \describe{
#'  \item{`matTumor`}{Matrix of log R ratio values *features x cells* for the
#'  tumor/sample cells (`matrix`).}
#'  \item{`matRef`}{Matrix of log R ratio values *features x cells* for the
#'  reference cells (`matrix`).}
#'  \item{`params`}{List of parameters set for the `genesPerWindow`, `refMeans`,
#'  `refMeanReads` and `thresh_capping` arguments (`list`).}
#'  \item{`coord`}{Data frame of coordinates for windows of peaks and associated
#'  data along the different steps (`data.frame`).
#'  Columns :
#'  - `CHR`, `start`, `end`, `id`: coordinates and name of genes.
#'  - `sumReads.tum/ref`: sum of read counts for all cells in tumor or reference
#'  cells.
#'  - `meanReads.tum/ref`: mean of read counts per cells for tumor or reference
#'  cells.
#'  - `sdReads.tum/ref`: standard deviation of read counts per cells for tumor
#'  or reference cells.
#'  - `keep`: logical, `TRUE` for genes to keep after filtering based on
#'  reference coverage (depends on `refReads` and `refMeanReads` arguments).
#'  - `meanReads/sdReads.norm.tum/ref`: mean/sd of normalized counts per million
#'  for tumor/reference cells.
#'  - `meanLRR/sdReads.raw.tum/ref`: mean/sd of raw log R ratio (LRR) for
#'  tumor/reference cells.
#'  - `meanLRR/sdLRR.cap.tum/ref`: mean/sd of capped log R ratio (LRR) for
#'  tumor/reference cells (depends on `thresh_capping` argument).
#'  - `meanLRR/sdLRR.smoo.tum/ref`: mean/sd of smoothed log R ratio (LRR) for
#'  tumor/reference cells (means of moving windows defined by the
#'  `genesPerWindow` argument).
#'  - `meanLRR/sdLRR.cent.tum/ref`: mean/sd of centered log R ratio (LRR) for
#'  tumor/reference cells.
#'  - `meanLRR/sdLRR.corr.tum/ref`: mean/sd of final log R ratio (LRR) corrected
#'  by reference variability for tumor/reference cells.}
#' }
#' If `all_steps` is set to `TRUE`, the previously described list can be found
#' for each step in a list element.
#'
#' @family computeLogRatio
#'
#' @import Matrix
#' @import dplyr
#' @importFrom GenomicRanges GRanges
#' @importFrom sparseMatrixStats rowSds
#' @importFrom caTools runmean
#'
#' @export
#'
#' @details
#' The raw count matrix is transformed into log R ratios through the following
#' steps:
#'  - Match genes in count matrix with coordinates
#'  - Filtering on coverage (`refReads` and `refMeanReads` arguments)
#'  - Normalization for sequencing depth
#'  - Log transformation and normalization by reference data: log R ratio
#'  - Capping the range of values (`thresh_capping` argument)
#'  - Smoothing on genes windows
#'  - Centering of cells
#'  - Correcting by reference variability
#'
#' @examples
#' obj_rna <- computeLogRatioRNA(
#'   matTumor = matCounts(muscadet)$RNA,
#'   matRef = matCounts(muscadet_ref)$RNA,
#'   genesCoord = coordFeatures(muscadet)$RNA,
#'   genome = slot(muscadet, "genome"),
#'   refReads = 20 # low value only for the example subsampled datasets
#' )
#' table(obj_rna$coord$keep)
#'
#' # With results form every step when `all_steps = TRUE`
#' obj_rna_all <- computeLogRatioRNA(
#'   matTumor = matCounts(muscadet)$RNA,
#'   matRef = matCounts(muscadet_ref)$RNA,
#'   genesCoord = coordFeatures(muscadet)$RNA,
#'   genome = slot(muscadet, "genome"),
#'   refReads = 20, # low value only for the example subsampled datasets
#'   all_steps = TRUE
#' )
#' names(obj_rna_all)
#' table(obj_rna_all$coord$keep)
#' nrow(obj_rna_all$step08$matTumor)
#'
computeLogRatioRNA <- function(matTumor,
                               matRef,
                               genesCoord,
                               genome = "hg38",
                               genesPerWindow = 101,
                               refReads = 100,
                               refMeanReads = 0.01,
                               thresh_capping = 3,
                               all_steps = FALSE,
                               quiet = FALSE) {
  params <- list(
    genesPerWindow = genesPerWindow,
    refReads = refReads,
    refMeanReads = refMeanReads,
    thresh_capping = thresh_capping
  )

  # Sort cells (column names of matrices)
  matTumor <- matTumor[, sort(colnames(matTumor))]
  matRef <- matRef[, sort(colnames(matRef))]

  ## Step 01: Match genes in count matrix with coordinates ---------------------
  if (quiet == FALSE) {
    message("Step 01 - Match genes in count matrix with coordinates")
  }

  # Initialize gene coordinates data frame
  genesDF <- as.data.frame(genesCoord)
  colnames(genesDF)[1] <- "CHR"
  genesDF$CHR <- factor(genesDF$CHR, levels = unique(genesDF$CHR))

  # Remove genes from matrix without coordinates in the input gene list
  matTumor <- matTumor[!is.na(match(rownames(matTumor), genesDF[, 4])), ]
  # Get same genes in reference data
  matTumor <- matTumor[rownames(matTumor) %in% rownames(matRef), ] # remove gene without data in reference
  matRef <- matRef[rownames(matTumor), ] # same genes in correct order in reference

  # Match genes data with rows of matrix in correct order
  genesDF <- genesDF[na.omit(match(rownames(matTumor), genesDF[, 4])), ]

  obj <- list()

  if (all_steps == T) {
    obj$step01 <- list(
      matTumor = matTumor,
      matRef = matRef,
      name = "Counts per gene"
    )
  }

  ## Step 02: Filtering on coverage --------------------------------------------
  if (quiet == FALSE) {
    message("Step 02 - Filtering genes: Minimum of ", refReads, " read(s) in reference cells and minimum of ", refMeanReads, " read(s) in average per reference cell")
  }

  # Add data to gene coordinates data frame
  genesDF$sumReads.tum <- Matrix::rowSums(matTumor)
  genesDF$meanReads.tum <- Matrix::rowMeans(matTumor)
  genesDF$sdReads.tum <- sparseMatrixStats::rowSds(matTumor)
  genesDF$sumReads.ref <- Matrix::rowSums(matRef)
  genesDF$meanReads.ref <- Matrix::rowMeans(matRef)
  genesDF$sdReads.ref <- sparseMatrixStats::rowSds(matRef)

  # Add filtering information 'keep' in gene coordinates data frame
  genesDF <- mutate(genesDF, keep = meanReads.ref >= refMeanReads & sumReads.ref >= refReads)

  stopifnot("No genes passing filters (`refReads` and `refMeanReads` arguments)" = any(genesDF$keep))

  # Filter matrices
  matTumor <- matTumor[genesDF$keep, ]
  matRef <- matRef[genesDF$keep, ]

  if (all_steps == T) {
    obj$step02 <- list(
      matTumor = matTumor,
      matRef = matRef,
      name = "Counts per gene (filtered)"
    )
  }

  ## Step 03: Normalization for sequencing depth -------------------------------
  if (quiet == FALSE) {
    message("Step 03 - Normalization for sequencing depth: Normalized counts per million")
  }

  # Normalization and transformation in counts per million
  matTumor <- Matrix(apply(matTumor, 2, function(x) {
    (x / sum(x, na.rm = T)) * 1e6
  }), sparse = TRUE)
  matRef <- Matrix(apply(matRef, 2, function(x) {
    (x / sum(x, na.rm = T)) * 1e6
  }), sparse = TRUE)

  if (all_steps == T) {
    obj$step03 <- list(
      matTumor = matTumor,
      matRef = matRef,
      name = "Normalized counts per million"
    )
  }

  genesDF$meanReads.norm.tum[genesDF$keep] <- Matrix::rowMeans(matTumor)
  genesDF$sdReads.norm.tum[genesDF$keep] <- sparseMatrixStats::rowSds(matTumor)
  genesDF$meanReads.norm.ref[genesDF$keep] <- Matrix::rowMeans(matRef)
  genesDF$sdReads.norm.ref[genesDF$keep] <- sparseMatrixStats::rowSds(matRef)


  ## Step 04: Log transformation and normalization by reference data: log R ratio ----------
  if (quiet == FALSE) {
    message("Step 04 - Log transformation and normalization by reference data: log R ratio")
  }

  # Log transformation
  matTumor <- Matrix(apply(matTumor, 2, function(x) {
    log2(1 + x)
  }), sparse = TRUE)
  matRef <- Matrix(apply(matRef, 2, function(x) {
    log2(1 + x)
  }), sparse = TRUE)

  # Ratio versus the mean of reference cells
  matTumor <- apply(matTumor, 2, "-", Matrix::rowMeans(matRef))
  matRef <- apply(matRef, 2, "-", Matrix::rowMeans(matRef))

  if (all_steps == T) {
    obj$step04 <- list(
      matTumor = matTumor,
      matRef = matRef,
      name = "Log R ratio"
    )
  }

  genesDF$meanLRR.raw.tum[genesDF$keep] <- Matrix::rowMeans(matTumor)
  genesDF$sdLRR.raw.tum[genesDF$keep] <- matrixStats::rowSds(matTumor)
  genesDF$meanLRR.raw.ref[genesDF$keep] <- Matrix::rowMeans(matRef)
  genesDF$sdLRR.raw.ref[genesDF$keep] <- matrixStats::rowSds(matRef)


  ## Step 05: Capping the range of values --------------------------------------
  if (quiet == FALSE) {
    message("Step 05 - Capping the range of values: threshold = ", thresh_capping)
  }

  # Max and min thresholds of log R ratios
  matTumor[matTumor > thresh_capping] <- thresh_capping
  matTumor[matTumor < (-1 * thresh_capping)] <- -1 * thresh_capping
  matRef[matRef > thresh_capping] <- thresh_capping
  matRef[matRef < (-1 * thresh_capping)] <- -1 * thresh_capping

  if (all_steps == T) {
    obj$step05 <- list(
      matTumor = matTumor,
      matRef = matRef,
      name = "Log R ratio capped"
    )
  }

  genesDF$meanLRR.cap.tum[genesDF$keep] <- Matrix::rowMeans(matTumor)
  genesDF$sdLRR.cap.tum[genesDF$keep] <- matrixStats::rowSds(matTumor)
  genesDF$meanLRR.cap.ref[genesDF$keep] <- Matrix::rowMeans(matRef)
  genesDF$sdLRR.cap.ref[genesDF$keep] <- matrixStats::rowSds(matRef)


  ## Step 06: Smoothing on genes windows ---------------------------------------
  if (quiet == FALSE) {
    message("Step 06 - Smoothing values on gene windows: ", genesPerWindow, " genes per window")
  }

  matTumor <- do.call(rbind, sapply(unique(genesDF[, 1]), function(i) {
    # Subset matrix per chromosomes
    mat <- matTumor[genesDF[genesDF[, 1] == i & genesDF$keep == TRUE, 4], ]
    if (nrow(mat) != 0) {
      row_names <- rownames(mat)
      col_names <- colnames(mat)
      # If the number of genes per window is higher than the number of genes on chromosome
      if (nrow(mat) < genesPerWindow) genesPerWindow <- nrow(mat)
      # Get moving window means
      mat_mean <- caTools::runmean(mat, k = genesPerWindow, align = "center")
      rownames(mat_mean) <- row_names
      colnames(mat_mean) <- col_names
      return(mat_mean)
    } else {
      return(NULL)
    }
  }))

  matRef <- do.call(rbind, sapply(unique(genesDF[, 1]), function(i) {
    # Subset matrix per chromosomes
    mat <- matRef[genesDF[genesDF[, 1] == i & genesDF$keep == TRUE, 4], ]
    if (nrow(mat) != 0) {
      row_names <- rownames(mat)
      col_names <- colnames(mat)
      # If the number of genes per window is higher than the number of genes on chromosome
      if (nrow(mat) < genesPerWindow) genesPerWindow <- nrow(mat)
      # Get moving window means
      mat_mean <- caTools::runmean(mat, k = genesPerWindow, align = "center")
      rownames(mat_mean) <- row_names
      colnames(mat_mean) <- col_names
      return(mat_mean)
    } else {
      return(NULL)
    }
  }))

  if (all_steps == T) {
    obj$step06 <- list(
      matTumor = matTumor,
      matRef = matRef,
      name = "Log R ratio smoothed"
    )
  }

  genesDF$meanLRR.smoo.tum[genesDF$keep] <- Matrix::rowMeans(matTumor)
  genesDF$sdLRR.smoo.tum[genesDF$keep] <- matrixStats::rowSds(matTumor)
  genesDF$meanLRR.smoo.ref[genesDF$keep] <- Matrix::rowMeans(matRef)
  genesDF$sdLRR.smoo.ref[genesDF$keep] <- matrixStats::rowSds(matRef)


  ## Step 07: Centering of cells ------------------------------
  if (quiet == FALSE) {
    message("Step 07 - Centering of cells")
  }

  # Center cells
  meansTumor <- apply(matTumor, 2, function(x) {
    mean(x, na.rm = TRUE)
  })
  matTumor <- t(apply(matTumor, 1, "-", meansTumor))

  meansRef <- apply(matRef, 2, function(x) {
    mean(x, na.rm = TRUE)
  })
  matRef <- t(apply(matRef, 1, "-", meansRef))

  if (all_steps == T) {
    obj$step07 <- list(
      matTumor = matTumor,
      matRef = matRef,
      name = "Log R ratio centered"
    )
  }

  genesDF$meanLRR.cent.tum[genesDF$keep] <- Matrix::rowMeans(matTumor)
  genesDF$sdLRR.cent.tum[genesDF$keep] <- matrixStats::rowSds(matTumor)
  genesDF$meanLRR.cent.ref[genesDF$keep] <- Matrix::rowMeans(matRef)
  genesDF$sdLRR.cent.ref[genesDF$keep] <- matrixStats::rowSds(matRef)


  ## Step 08: Correcting by reference variability ------------------------------
  if (quiet == FALSE) {
    message("Step 08 - Correcting by reference variability")
  }

  sd.ref <- apply(matRef, 1, function(x) sd(x, na.rm = TRUE))
  mean.ref <- apply(matRef, 1, function(x) mean(x, na.rm = TRUE))

  matTumor <- apply(matTumor, 2, function(x) {
    (x - mean.ref) / sd.ref
  })
  matRef <- apply(matRef, 2, function(x) {
    (x - mean.ref) / sd.ref
  })

  if (all_steps == T) {
    obj$step08 <- list(
      matTumor = matTumor,
      matRef = matRef,
      name = "Log R ratio corrected by ref variablity"
    )
  }

  genesDF$meanLRR.corr.tum[genesDF$keep] <- Matrix::rowMeans(matTumor)
  genesDF$sdLRR.corr.tum[genesDF$keep] <- matrixStats::rowSds(matTumor)
  genesDF$meanLRR.corr.ref[genesDF$keep] <- Matrix::rowMeans(matRef)
  genesDF$sdLRR.corr.ref[genesDF$keep] <- matrixStats::rowSds(matRef)

  if (all_steps == T) {
    obj$params <- params
    obj$coord <- genesDF
    return(obj)
  } else {
    obj_min <- list(
      matTumor = matTumor,
      matRef = matRef,
      params = params,
      coord = genesDF
    )
    return(obj_min)
  }
}
