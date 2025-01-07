#' Assign a cluster assignment to a `muscadet` object
#'
#' Add the user selected cluster assignments to cells in a
#' \code{\link{muscadet}} object. This function allows the user to choose the
#' cluster assignments they consider to fit the data and their requirements, or
#' cluster assignments based on other data and methods.
#'
#' @param x A \code{\link{muscadet}} object (`muscadet`).
#' @param k Integer specifying the number of clusters to choose from the
#'   muscadet object (`integer`). If `k` is provided, the clustering result
#'   stored in the object for that specific `k` is used
#'   (`muscadet_obj$clustering$clusters[[as.character(k)]]`).
#' @param clusters A custom named vector of cluster assignments (e.g., to
#'   regroup some clusters) (`vector`). The vector must have the same length as
#'   the number of cells in the muscadet object and names matching the cell
#'   names.
#' @param mapping Optional named vector specifying how to remap cluster values
#'   (`vector`). The names of the vector correspond to the original cluster
#'   values, and the values are the remapped cluster values. For example, `c("1"
#'   = 1, "2" = 1, "3" = 2, "4" = 3)` would merge clusters 1 and 2 into 1,
#'   cluster 3 into 2, and cluster 4 into 3.
#'
#' @details
#' - The clusters can be taken directly from the `muscadet` object clustering
#' results with setting the `k` argument to the chosen partition (e.g.
#' `muscadet_obj$clustering$clusters[["2"]]` for k=`2`).
#' - A custom vector of cluster assignments
#' can be attributed using the `clusters` argument.
#' - Either way, the clusters assignments can be rearranged using the `mapping`
#' argument.
#'
#' @return A \code{\link{muscadet}} object updated with the user chosen cluster
#' assignments in `muscadet_obj$cnacalling$clusters`.
#'
#' @include objects.R
#' @importFrom SeuratObject Cells
#' @export
#'
#' @examples
#' # Load muscadet object
#' data(muscadet_obj)
#'
#' # Select clustering result for k = 3
#' muscadet_obj <- assignClusters(muscadet_obj, k = 3)
#' table(muscadet_obj$cnacalling$clusters)
#'
#' # Assign custom clusters
#' cell_names <- Reduce(union, SeuratObject::Cells(muscadet_obj))
#' n1 <- sample(1:length(cell_names), 1)
#' n2 <- length(cell_names) - n1
#' custom_clusters <- c(rep.int(1, n1), rep.int(2, n2))
#' names(custom_clusters) <- cell_names
#' muscadet_obj <- assignClusters(muscadet_obj, clusters = custom_clusters)
#' table(muscadet_obj$cnacalling$clusters)
#'
#' # Assign clusters with remapping
#' # example to remap from k=4 to k=3 by merging clusters 1 and 2
#' clusters <- muscadet_obj$clustering$clusters[["4"]]
#' mapping <- c("1" = 1, "2" = 1, "3" = 2, "4" = 3)
#'
#' muscadet_obj <- assignClusters(muscadet_obj, clusters = clusters, mapping = mapping)
#' table(muscadet_obj$cnacalling$clusters)
#' # check original and remapped clusters
#' table(clusters, muscadet_obj$cnacalling$clusters)
#'
#' muscadet_obj <- assignClusters(muscadet_obj, k=4, mapping = mapping)
#' table(muscadet_obj$cnacalling$clusters)
#' # check original and remapped clusters
#' table(clusters, muscadet_obj$cnacalling$clusters)
#'
#' # Visualize clusters on heatmap
#' ht <- heatmapMuscadet(
#'   muscadet_obj,
#'   k = 3,
#'   filename = file.path("heatmap_muscadet_k3.png"),
#'   title = "Example sample | k=3"
#' )
#' ht <- heatmapMuscadet(
#'   muscadet_obj,
#'   clusters = muscadet_obj$cnacalling$clusters,
#'   filename = file.path("heatmap_muscadet_customk3.png"),
#'   title = "Example sample | rearranged clusters k=3"
#' )
#'
assignClusters <- function(x, k = NULL, clusters = NULL, mapping = NULL) {

    # Validate input: x must be a muscadet object
    stopifnot("Input object 'x' must be of class 'muscadet'." = inherits(x, "muscadet"))

    # Validate that either k or clusters is provided, but not both
    stopifnot(
        "Either 'k' or 'clusters' must be provided, but not both." = xor(!is.null(k), !is.null(clusters))
    )

    # If k is provided, validate that clustering has been performed
    if (!is.null(k)) {
        stopifnot(
            "Clustering results are not available in the muscadet object. Perform clustering first using 'clusterMuscadet()'." = !is.null(x@clustering),
            "Clustering results for the chosen k are not available. Use 'clusterMuscadet()' with the 'k_range' argument containing k" = as.character(k) %in% names(x@clustering$clusters)
        )
    }

    # If clusters is provided, validate its format
    if (!is.null(clusters)) {
        stopifnot(
            "'clusters' must be a named vector." = is.vector(clusters) && !is.null(names(clusters)),
            "Length of 'clusters' must match the total number of cells in the muscadet object." = length(clusters) == length(Reduce(union, Cells(x))),
            "Names of 'clusters' must match the cell names in the muscadet object." = all(names(clusters) %in% Reduce(union, Cells(x)))
        )
    }

    # Apply mapping if provided
    if (!is.null(mapping)) {
        stopifnot(
            "Mapping must be a named vector." = is.vector(mapping) && !is.null(names(mapping)),
            "All cluster values must have corresponding mappings." = all(as.character(unique(clusters)) %in% names(mapping))
        )
        if (!is.null(k)) {
            clusters <- x@clustering$clusters[[as.character(k)]]
        }
        remapped_clusters <- mapping[as.character(clusters)]
        names(remapped_clusters) <- names(clusters)
        x@cnacalling[["clusters"]] <- remapped_clusters
    } else {
        # Assign clusters without remapping
        if (!is.null(clusters)) {
            x@cnacalling[["clusters"]] <- clusters
        } else if (!is.null(k)) {
            x@cnacalling[["clusters"]] <- x@clustering$clusters[[as.character(k)]]
        }
    }

    return(x)
}


#' Add allele counts to a `muscadet` object
#'
#' This function adds allele counts data to the `omics` of a
#' \code{\link{muscadet}} object. The data frames in the `allele_counts` list
#' are assigned to the `allelic` slots of `omics`.
#'
#' @param x A \code{\link{muscadet}} object.
#' @param allele_counts A list of data frames where each data frame contains
#'   allele counts for a specific omic (`list`). The list must have the same
#'   length and order as the number of omics in the `x` object. Each data frames
#'   must contain the following columns : `cell`, `id`, `CHROM`, `POS`, `REF`,
#'   `ALT`, `RD`, `AD`, `DP`, `GT`. See [allele_counts] for details.
#'
#' @return
#' A modified \code{\link{muscadet}} object with updated allele counts in the
#' `allelic` slot of each `muscomic` object in the `omics` slot.
#'
#' @note
#' As the allele counts are not used for computing log R
#' ratios and cell clustering, they are not mandatory at the creation of
#' `muscomic` and `muscadet` objects. The allele counts data can be added to
#' objects later with this `addAlleleCounts` function, before using the
#' [muscadet::mergeCounts()] function.
#'
#' This function is also useful to add allele counts for sample-specific SNP
#' positions to a common `muscadet` object, for example for the reference
#' `muscadet` object: a common `muscadet` object with reference coverage data
#' can be stored as a unique object to use as reference for computing log R
#' ratios for all samples, and then it can updated with allele counts at sample
#' specific SNP positions (e.g. found by bulk sequencing) before Copy Number
#' Alterations (CNAs) calling.
#'
#' @seealso [muscadet::CreateMuscomicObject()], [muscadet::mergeCounts()]
#'
#' @importFrom stringr str_remove
#' @importFrom gtools mixedsort
#'
#' @export
#'
#' @examples
#' # Load muscadet objects
#' data(muscadet_obj)
#' data(muscadet_obj_ref)
#'
#' # Add allele counts data frames to muscadet objects
#' muscadet_obj <- addAlleleCounts(
#'     muscadet_obj,
#'     allele_counts = list(allele_counts_atac_tumor, allele_counts_rna_tumor))
#' muscadet_obj_ref <- addAlleleCounts(
#'     muscadet_obj_ref,
#'     allele_counts = list(allele_counts_atac_ref, allele_counts_rna_ref))
#'
addAlleleCounts <- function(x, allele_counts) {

    # Check that x is a muscadet object
    stopifnot("The argument 'x' must be a muscadet object." = inherits(x, "muscadet"))

    # Check that allele_counts is a list of dataframes
    stopifnot("The argument 'allele_counts' must be a list." = inherits(allele_counts, "list"))
    stopifnot("The argument 'allele_counts' must be a list of data frames." = all(unlist(
        lapply(allele_counts, function(x)
            inherits(x, "data.frame"))
    )))

    # Ensure allele_counts list matches the number of omics in x
    if (length(allele_counts) != length(x@omics)) {
        stop("The 'allele_counts' list must have the same length as the number of omics in 'x'.")
    }

    # Assign names to allele_counts based on omics names in x
    names(allele_counts) <- names(x@omics)
    # Extract types from omics in x
    type.omic <- unlist(lapply(x@omics, function(x) x$type))

    # Format allele count data frames
    allele_counts <- lapply(allele_counts, function(allele_df) {
        # Remove "chr" if necessary
        allele_df$CHROM <- stringr::str_remove(allele_df$CHROM, "chr")
        # Ordered chromosomes
        allele_df$CHROM <- ordered(allele_df$CHROM, levels = gtools::mixedsort(unique(allele_df$CHROM)))
        # Reorder data
        allele_df <- allele_df[order(allele_df[, "POS"]), ]
        allele_df <- allele_df[order(allele_df[, "CHROM"]), ]
        return(allele_df)

    })

    # Add allele counts to each omic in the muscadet object
    for (i in names(x@omics)) {
        slot(x@omics[[i]], "allelic") <- list(table.counts = data.frame(omic = type.omic[[i]], allele_counts[[i]]))
    }

    return(x)
}


#' Merge counts for `muscadet` objects
#'
#' This function combines allelic and coverage counts from all omics per cluster
#' for both sample and reference. The resulting merged data is
#' stored in the `cnacalling` slot of the sample `muscadet` object.
#'
#' @param x A \code{\link{muscadet}} object containing sample data (`muscadet`).
#'   This object must include clustering assignments in the
#'   `cnacalling$clusters` slot.
#' @param reference A \code{\link{muscadet}} object containing reference data
#'   (`muscadet`).
#'
#' @return
#' A modified \code{\link{muscadet}} object corresponding to the `x` muscadet object,
#' with updated `cnacalling` slot containing:
#' \itemize{
#'   \item \code{allelic.counts}: Processed allelic counts on SNP positions, for all omics
#'   \item \code{coverage.counts}: Processed coverage counts merged with the reference.
#'   \item \code{combined.counts}: Combined data for allelic and coverage counts.
#' }
#' Abbreviations:
#' - RD = Reference allele read depth
#' - AD = Alternative allele read depth
#' - DP = Total read depth
#' - TUM = tumor sample
#' - NOR = normal reference
#' - omic = omic specific (`omic` column)
#' - all = for all omics
#'
#' @import dplyr
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' # Load muscadet objects
#' data(muscadet_obj)
#' data(muscadet_obj_ref)
#' # Merge counts from all omics from both sample and reference
#' muscadet_obj <- mergeCounts(muscadet_obj, muscadet_obj_ref)
#'
mergeCounts <- function(x, reference) {

    # Validate input: x and reference must be muscadet objects
    stopifnot("Input object 'x' must be of class 'muscadet'." = inherits(x, "muscadet"))
    stopifnot("Input object 'reference' must be of class 'muscadet'." = inherits(reference, "muscadet"))

    # Ensure clustering data is available in the sample muscadet object
    stopifnot(
        "Cluster assignments not found in the 'x' muscadet object. Use 'assignClusters()' to add them."
        = !is.null(x@cnacalling$clusters)
    )

    # Ensure allele counts data is available in the sample and reference muscadet object
    stopifnot(
        "Input object 'x' must contain allele counts data in the 'allelic' slot of omics. Use 'addAlleleCounts()' to add them."
        = all(unlist(lapply(x@omics, function(omic) {
            !is.null(omic@allelic$table.counts)
        })))
    )
    stopifnot(
        "Input object 'reference' must contain allele counts data in the 'allelic' slot of omics. Use 'addAlleleCounts()' to add them."
        = all(unlist(lapply(reference@omics, function(omic) {
            !is.null(omic@allelic$table.counts)
        })))
    )

    # Allelic Data Processing --------------------------------------------------
    # Extract allelic counts for all omics in both sample and reference
    allele_sample <- Reduce(rbind, lapply(x@omics, function(omic) {
        omic@allelic$table.counts
    }))
    allele_ref <- Reduce(rbind, lapply(reference@omics, function(omic) {
        omic@allelic$table.counts
    }))

    # Filter for cells with valid cluster assignments
    allele_sample <- allele_sample[allele_sample$cell %in% names(x@cnacalling$clusters),]
    # Add cluster assignments to table
    allele_sample$cluster <- x@cnacalling$clusters[match(allele_sample$cell, names(x@cnacalling$clusters))]

    # Calculate per-omic and total counts for the sample
    allele_sample <- allele_sample %>%
        dplyr::group_by(.data$cluster, .data$id, .data$omic) %>%
        dplyr::mutate(
            RD.omic = sum(.data$RD),
            AD.omic = sum(.data$AD),
            DP.omic = sum(.data$DP),
            AF.omic = round(.data$AD.omic / .data$DP.omic, 2)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$cluster, .data$id) %>%
        dplyr::mutate(
            RD.all = sum(.data$RD, na.rm = T),
            AD.all = sum(.data$AD, na.rm = T),
            DP.all = sum(.data$DP, na.rm = T),
            AF.all = round(.data$AD.all / .data$DP.all, 2)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(!c("cell", "RD", "AD", "DP"))
    allele_sample <- subset(allele_sample, !duplicated(allele_sample))

    # Calculate per-omic and total counts for the reference
    allele_ref <- allele_ref %>%
        dplyr::group_by(.data$id, .data$omic) %>%
        dplyr::mutate(
            RD.omic = sum(.data$RD),
            AD.omic = sum(.data$AD),
            DP.omic = sum(.data$DP),
            AF.omic = round(.data$AD.omic / .data$DP.omic, 2)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$id) %>%
        dplyr::mutate(
            RD.all = sum(.data$RD, na.rm = T),
            AD.all = sum(.data$AD, na.rm = T),
            DP.all = sum(.data$DP, na.rm = T),
            AF.all = round(.data$AD.all / .data$DP.all, 2)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(!c("cell", "RD", "AD", "DP"))
    allele_ref <- subset(allele_ref, !duplicated(allele_ref))

    # Merge sample and reference allelic data
    snp <- dplyr::left_join(
        allele_sample,
        allele_ref,
        by = c("CHROM", "POS", "REF", "ALT", "id", "GT", "omic"),
        suffix = c(".TUM", ".NOR")
    ) %>%
        dplyr::arrange(.data$CHROM, .data$POS, .data$cluster) %>%
        dplyr::mutate(signal = "allelic")

    x@cnacalling[["allelic.counts"]] <- as.data.frame(snp)

    # Coverage Data Processing -------------------------------------------------
    # Extract coverage counts for all omics in both sample and reference
    coverage_sample <- Reduce(rbind, lapply(x@omics, function(omic) {
        omic@coverage$table.counts
    }))
    coverage_ref <- Reduce(rbind, lapply(reference@omics, function(omic) {
        omic@coverage$table.counts
    }))

    # Filter for cells with valid cluster assignments
    coverage_sample <- coverage_sample[coverage_sample$cell %in% names(x@cnacalling$clusters),]
    # Add cluster assignments to table
    coverage_sample$cluster <- x@cnacalling$clusters[match(coverage_sample$cell, names(x@cnacalling$clusters))]

    # Calculate total counts for the sample
    coverage_sample <- coverage_sample %>%
        dplyr::group_by(.data$cluster, .data$id) %>%
        dplyr::mutate(DP = sum(.data$DP)) %>%
        dplyr::ungroup() %>%
        dplyr::select(!c("cell"))
    coverage_sample <- subset(coverage_sample, !duplicated(coverage_sample))

    # Calculate total counts for the reference
    coverage_ref <- coverage_ref %>%
        dplyr::group_by(.data$id) %>%
        dplyr::mutate(DP = sum(.data$DP)) %>%
        dplyr::ungroup() %>%
        dplyr::select(!c("cell"))
    coverage_ref <- subset(coverage_ref, !duplicated(coverage_ref))

    # Merge sample and reference coverage data
    cov <- dplyr::left_join(
        coverage_sample,
        coverage_ref,
        by = c("CHROM", "POS", "id", "omic"),
        suffix = c(".TUM", ".NOR")
    ) %>%
        dplyr::arrange(.data$CHROM, .data$POS, .data$cluster) %>%
        dplyr::select("omic", "id", "CHROM", "POS", "DP.NOR", "DP.TUM", "cluster") %>%
        dplyr::mutate(signal = "coverage")

    x@cnacalling[["coverage.counts"]] <- as.data.frame(cov)


    # Combine Allelic and Coverage Data ----------------------------------------

    # Format data
    snp <- snp[, c("CHROM", "POS", "DP.all.NOR", "RD.all.NOR", "DP.all.TUM", "RD.all.TUM", "cluster", "signal")]
    snp <- unique(snp)
    cov <- cov[, c("CHROM", "POS", "DP.NOR", "DP.NOR", "DP.TUM", "DP.TUM", "cluster", "signal")]
    cov <- unique(cov)

    colnames(snp) <- c("Chromosome", "Position", "NOR.DP", "NOR.RD", "TUM.DP", "TUM.RD", "cluster", "signal")
    colnames(cov) <- colnames(snp)

    # snp$NOR.RD = round(snp$NOR.RD / 2, 0)

    # Make sure the levels match for binding
    levels(snp$Chromosome) <- union(levels(snp$Chromosome), levels(cov$Chromosome))
    levels(cov$Chromosome) <- union(levels(snp$Chromosome), levels(cov$Chromosome))

    # Combine both data
    combined <- rbind(snp, cov) %>%
        arrange(.data$Chromosome, .data$Position)

    x@cnacalling[["combined.counts"]] <- as.data.frame(combined)

    return(x)
}
