#' Heatmap Plot Function for muscadet Object
#'
#' This function generates a heatmap to visualize log R ratio (LRR) data
#' contained in muscadet objects. One heatmap is generated per omic, rows are
#' cells and columns are chromosomes, for muscadet object containing multiple
#' omics, the heatmaps are plotted horizontally aligned. The cells can be
#' clustered for a specific k number of cluster following the clustering step of
#' muscadet object, or custom cluster assignments can be used. Additionally, LRR
#' values from bulk sequencing data can be plotted as an annotation under the
#' heatmaps.
#'
#' @param x A muscadet object containing LRR data for all omics (with
#'   [muscadet::computeLogRatio()]) and clustering data (with
#'   [muscadet::clusterMuscadet()]) (`muscadet`).
#'
#' @param filename Character string specifying the file path to save the heatmap
#'   image in the PNG format (`character` string).
#'
#' @param k (Optional) Integer specifying the cluster number to plot
#'   (`integer`). It should be within the range of k used for clustering (with
#'   [muscadet::clusterMuscadet()]). Should be provided if `clusters` is `NULL`.
#'
#' @param clusters (Optional) A custom named vector of cluster assignments
#'   (`integer` named vector). Names must corresponds the cell names of the
#'   `muscadet` object: to the totality of cells if `show_missing = TRUE`, or to
#'   only common cells (cells with data in all omics) if `show_missing = FALSE`.
#'   Should be provided if `k` is `NULL`.
#'
#' @param title Character string for the title of the heatmap (`character` string).
#'
#' @param add_bulk_lrr Logical value indicating whether to add bulk LRR
#' data as annotations if present in the muscadet object (`logical`). Default is `TRUE`.
#'
#' @param show_missing Logical value indicating whether to show missing cells
#' (cells with missing data in at least one omic) in the heatmaps (`logical`). Default is `TRUE`.
#'
#' @param white_scale Numeric vector of length 2 with values between 0 and 1,
#'   defining the white color boundaries (`numeric` vector). This parameter
#'   specifies the quantiles of the LRR reference data that will define
#'   the boundaries for the white color in the heatmap. LRR values falling
#'   within this range are considered close enough to the majority of the LRR
#'   reference data, indicating no significant gain or loss of coverage, thereby
#'   are represented as white on the heatmap.
#'
#'
#' @param colors Vector of colors for the cluster annotation. Default is `NULL`, which uses predefined colors.
#'
#' @param quiet Logical value indicating whether to suppress messages. Default is `FALSE`.
#'
#' @return A heatmap plot object saved as a png image using the provided file name.
#'
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#' @importFrom methods slot
#' @importFrom stats median
#' @importFrom grDevices pdf png palette dev.off
#' @importFrom grid gpar unit grid.rect grid.text
#'
#' @export
#'
#' @examples
#' # Load example muscadet object
#' data(muscadet_obj)
#'
#' # Perform clustering (if not already done)
#' muscadet_obj <- clusterMuscadet(muscadet_obj, k_range = 2:5)
#'
#' # Custom title with clustering parameters
#' title <- paste("Example sample |",
#'      slot(muscadet_obj, "clustering")[["params"]][["dist_method"]],
#'      slot(muscadet_obj, "clustering")[["params"]][["hclust_method"]], "|",
#'      "k=3", "|",
#'      "Weights of omics:",
#'      paste(slot(muscadet_obj, "clustering")[["params"]][["weights"]], collapse = ", "))
#'
#' # Generate heatmap
#' heatmapMuscadet(
#'     muscadet_obj,
#'     k = 3,
#'     filename = file.path(getwd(), "heatmap_muscadet_k3.png"),
#'     title = title
#' )
#'
#' heatmapMuscadet(
#'     muscadet_obj,
#'     k = 3,
#'     filename = file.path(getwd(), "heatmap_muscadet_k3_commoncells.png"),
#'     title = title,
#'     show_missing = FALSE
#' )
#'
heatmapMuscadet <- function(x, filename, k = NULL, clusters = NULL, title = "",
                            add_bulk_lrr = TRUE, show_missing = TRUE, white_scale = c(0.3, 0.7),
                            colors = NULL, quiet = FALSE) {

    # Validate muscadet object
    stopifnot(
        "The muscadet object does not contain clustering data (use clusterMuscadet() to perform clustering of log R ratio data)." =
            !is.null(slot(x, "clustering"))
    )

    # Validate k and clusters
    stopifnot(
        "Both k and clusters cannot be NULL." = !(is.null(k) && is.null(clusters))
    )

    # Validate k in clustering slot
    if (!is.null(k)) {
        stopifnot(
            "The muscadet object must contain the clustering results for the k provided." = as.character(k) %in% names(slot(x, "clustering")[["clusters"]])
        )
    }

    # Validate white_scale
    stopifnot(
        "white_scale must be a numeric vector of length 2." = length(white_scale) == 2
    )
    stopifnot(
        "white_scale values must be between 0 and 1." = all(white_scale >= 0 & white_scale <= 1)
    )
    white_scale <- round(white_scale, 2)
    stopifnot(
        "The two elements of white_scale must not be equal." = white_scale[1] != white_scale[2]
    )
    if (white_scale[1] > white_scale[2]) {
        message("The first element of white_scale is greater than the second. Reversing the values.")
        white_scale <- sort(white_scale)
    }

    # Set default color palette for clusters if not provided
    if (is.null(colors)) {
        colors <- c("#FABC2A", "#7FC97F", "#EE6C4D", "#39ADBD", "#BEAED4",
                              "#FEE672", "#F76F8E", "#487BEA", "#B67BE6", "#F38D68",
                              "#7FD8BE", "#F2AFC0")
    }
    # Set palette
    palette(colors)

    # Check if output directory exists
    stopifnot("The directory doesn't exist" = file.exists(dirname(filename)))

    # Get common and all cells
    common_cells <- sort(Reduce(intersect, lapply(muscadet::matLogRatio(x), colnames)))
    all_cells <- sort(Reduce(union, lapply(muscadet::matLogRatio(x), colnames)))

    if (quiet == FALSE) {
        # Print information messages
        message("---- Heatmap Parameters ----")
        message("Omics in the muscadet object: ", paste(sapply(slot(x, "omics"), function(n) slot(n, "label.omic")), collapse = ", "))
        omic_dims <- sapply(slot(x, "omics"), function(m) {
            paste0(m@label.omic, ": ", nrow(muscadet::matLogRatio(m)), " cells x ", ncol(muscadet::matLogRatio(m)), " features")
        })
        message("Omics log R ratio data dimensions:\n  ", paste(omic_dims, collapse = "\n  "))

        if (!is.null(k))
            message("Number of clusters (k): ", k)
        if (!is.null(clusters))
            message("Custom clusters provided: ", length(unique(clusters)), " clusters.")
        message("Number of cells: ", length(all_cells), " total (", length(common_cells), " common across all omics).")

        message("Show missing cells: ", show_missing)
        message("Bulk LRR annotations: ", ifelse(add_bulk_lrr, slot(x, "bulk.data")[["label"]], add_bulk_lrr))
        message("White scale quantiles: ", paste(white_scale, collapse = " - "))
        message("Output file: ", filename)
    }

    # Create list of heatmap objects
    list_ht <- lapply(slot(x, "omics"), function(muscomic) {

        pdf(file = NULL)
        ht_opt(message = F)

        # Get chromosome info for features
        coord <- muscomic@coverage$coord.features
        chrom <- factor(coord[coord$keep, "CHR"], levels = unique(coord[coord$keep, "CHR"]))

        # Define color breaks for heatmap using white_scale argument
        col_breaks <- c(-5, muscomic@coverage$ref.log.ratio.perc[as.character(white_scale)], 5)

        # Extract log R ratio matrix
        if (show_missing == TRUE) {
            mat <- t(muscadet::matLogRatio(muscomic))
            cells.diff <- setdiff(all_cells, rownames(mat)) # identify missing cells
            if (length(cells.diff) > 0) {
                mat.na <- matrix(data = NA, nrow = length(cells.diff), ncol = ncol(mat),
                                 dimnames = list(cells.diff, colnames(mat)))
                mat <- rbind(mat, mat.na)[all_cells, ] # add empty rows to matrix
            } else {
                mat <- mat[all_cells, ]
            }
        }
        if (show_missing == FALSE) {
            mat <- t(muscadet::matLogRatio(muscomic))
            mat <- mat[common_cells, ]
        }

        # Create heatmap
        ht <- ComplexHeatmap::Heatmap(
            mat,
            name = muscomic@label.omic,
            heatmap_legend_param = list(title = muscomic@label.omic),
            show_column_names = F,
            show_row_names = F,
            cluster_columns = F,
            column_split = chrom,
            column_title = paste(muscomic@label.omic, "coverage on", muscomic@coverage$label.features),
            row_title_gp = gpar(fontsize = 10),
            column_title_gp = gpar(fontsize = 10),
            border_gp = gpar(col = "black", lwd = 1),
            heatmap_height = unit(12, "cm"),
            heatmap_width = unit(18, "cm"),
            col = circlize::colorRamp2(col_breaks, c("#00008E", "white", "white", "#630000")),
            row_title_rot = 0,
            raster_device = "png",
            raster_quality = 3
        )

        # Add empty chromosome annotation
        emp <- ComplexHeatmap::anno_empty(border = FALSE, height = unit(2, "mm"))
        ht <- ComplexHeatmap::attach_annotation(ht, columnAnnotation(
            emp_annot = emp,
            name = paste0("chrom_", muscomic@label.omic)
        ), side = "top")
        # Rename chromosome annotation
        ht@top_annotation@anno_list[["emp_annot"]]@name <- paste0("chrom_", muscomic@label.omic)
        ht@top_annotation@anno_list[["emp_annot"]]@label <- paste0("chrom_", muscomic@label.omic)
        ht@top_annotation@anno_list[["emp_annot"]]@name_param[["label"]] <- paste0("chrom_", muscomic@label.omic)
        names(ht@top_annotation@anno_list)[1] <- paste0("chrom_", muscomic@label.omic)

        # Add bulk LRR data as annotation
        if (add_bulk_lrr == TRUE) {
            # Retrieve bulk lrr values on features
            bulk_df <- muscadet::getLogRatioBulk(muscomic, x@bulk.data$log.ratio)
            # Define color scale
            bulk_col <- list(circlize::colorRamp2(c(min(bulk_df$bulk.lrr), median(bulk_df$bulk.lrr), max(bulk_df$bulk.lrr)),
                                                  c("#00008E", "white", "#630000")))
            names(bulk_col) = paste0("bulk_", muscomic@label.omic)
            # Create data frame for annotation
            bulk_df <- as.data.frame(bulk_df$bulk.lrr)
            colnames(bulk_df) <- paste0("bulk_", muscomic@label.omic)
            # Create bulk LRR annotation
            annot_bulk <- ComplexHeatmap::HeatmapAnnotation(
                df = bulk_df,
                col = bulk_col,
                annotation_label = x@bulk.data$label,
                annotation_name_side = "left",
                annotation_name_gp = gpar(cex = 0.7)
            )
            # Attach annotation to heatmap
            ht <- ComplexHeatmap::attach_annotation(ht, annot_bulk, side = "bottom")
        }

        dev.off()
        return(ht)
    })


    # Create empty annotation to add number of cells (ncells)
    ha <- ComplexHeatmap::rowAnnotation(ncells = anno_empty(border = FALSE, width = unit(12, "mm")))
    # Combine heatmaps and ncells empty annotation
    ht_list <- ha + Reduce("+", list_ht)

    # Create layout for the list of heatmaps
    pdf(file = NULL)
    if (show_missing == TRUE) {
        if (is.null(clusters)) {
            # 1. Cluster assignments from the muscadet object clustering with all cells
            clusters <- x@clustering$clusters[[as.character(k)]]
        }
        n_cells <- table(clusters)

        ht_all <- ComplexHeatmap::draw(ht_list, column_title = title, ht_gap = unit(1, "cm"),
                             row_split = as.factor(clusters), row_order = names(clusters),
                             cluster_rows = F, merge_legend = TRUE)


    } else if (show_missing == FALSE) {
        if (is.null(clusters)) {
            # 2. Clustering hclust from the muscadet object clustering with commom cells
            hc <- x@clustering$hclust # hclust object to print the dendrogram on the heatmap
            n_cells <- table(stats::cutree(hc, as.integer(k)))

            ht_all <- ComplexHeatmap::draw(ht_list, column_title = title, ht_gap = unit(1, "cm"),
                                 cluster_rows = hc, row_split = as.integer(k), merge_legend = TRUE)
        } else {
            # 3. Custom cluster assignments vector
            ht_all <- ComplexHeatmap::draw(ht_list, column_title = title, ht_gap = unit(1, "cm"),
                                 row_split = as.factor(clusters), row_order = names(clusters),
                                 cluster_rows = F, merge_legend = TRUE)
        }
    }
    dev.off()

    # Save complete plot of heatmaps as PNG
    png(filename = filename,
        width = ht_all@ht_list_param[["width"]],
        height = ht_all@ht_list_param[["height"]],
        units = "mm", res = 300)

    # Print plot
    print(ht_all)

    # Add annotation: number of cells per cluster
    for(i in 1:length(n_cells)) {
        ComplexHeatmap::decorate_annotation("ncells", slice = i, envir = environment(), {
            grid.rect(x = 1, width = unit(2, "mm"), gp = gpar(fill = colors[i], col = NA), just = "right")
            grid.text(n_cells[i], x = 0.7, just = "right", gp = gpar(cex = 0.75))
        })
    }

    # Add annotation: chromosome numbers
    for (ht_name in names(ht_all@ht_list)[-1]) {
        chrom_names <- unique(names(ht_all@ht_list[[ht_name]]@column_order_list))
        for (i in 1:length(chrom_names)) {
            ComplexHeatmap::decorate_annotation(paste0("chrom_", ht_name), slice = i, envir = environment(), {
                grid.text(chrom_names[i], just = "center", gp = gpar(fontsize=8.5))
            })
        }
    }

    dev.off()
}




