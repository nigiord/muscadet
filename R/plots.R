#' Heatmap plot for `muscadet` object
#'
#' This function generates a heatmap to visualize log R ratio (LRR) data
#' contained in \code{\link{muscadet}} objects. One heatmap is generated per
#' omic, rows are cells and columns are chromosomes, for `muscadet` object
#' containing multiple omics, the heatmaps are plotted horizontally aligned. The
#' cells can be clustered for a specific k number of cluster following the
#' clustering step of `muscadet` object, or custom cluster assignments can be
#' used. Additionally, LRR values from bulk sequencing data can be plotted as an
#' annotation under the heatmaps.
#'
#' @param x A \code{\link{muscadet}} object containing LRR data for all omics (with
#'   [muscadet::computeLogRatio()]) and clustering data (with
#'   [muscadet::clusterMuscadet()]) (`muscadet`).
#'
#' @param filename (Optional) Character string specifying the file path to save
#'   the heatmap image in the PNG format (`character` string).
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
#' @param title Character string for the title of the heatmap (`character`
#'   string). Default is an empty character string.
#'
#' @param add_bulk_lrr Logical value indicating whether to add bulk LRR data as
#'   annotations if present in the muscadet object (`logical`). By default, it
#'   is set to `NULL`, which will automatically include bulk LRR data if it is
#'   available in the muscadet object.
#'
#' @param show_missing Logical value indicating whether to show missing cells
#'   (cells with missing data in at least one omic) in the heatmaps (`logical`).
#'   Default is `TRUE`.
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
#' @param colors Vector of colors for the cluster annotation (`character`
#'   vector). Default is `NULL`, which uses predefined colors.
#'
#' @param quiet Logical value indicating whether to suppress messages. Default
#'   is `FALSE`.
#'
#' @return
#' A list containing:
#' - `plot`: A \code{\link{gTree}} object created with [grid::grid.grab()] (\code{\link{gTree}}).
#' - `width`: Width of the heatmap plot in mm (\code{\link{unit}}).
#' - `height`: Height of the heatmap plot in mm (\code{\link{unit}}).
#'
#' If the `filename` argument is provided, the heatmap is directly saved as a
#' PNG image at the provided path.
#'
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#' @importFrom methods slot
#' @importFrom stats median
#' @importFrom grDevices pdf png palette dev.off
#' @importFrom grid gpar unit grid.rect grid.text grid.grab
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
#'     filename = file.path("heatmap_muscadet_k3.png"),
#'     title = title
#' )
#'
#' heatmapMuscadet(
#'     muscadet_obj,
#'     k = 3,
#'     filename = file.path("heatmap_muscadet_k3_commoncells.png"),
#'     title = title,
#'     show_missing = FALSE
#' )
#'
#' # Loop over k
#' for (k in names(muscadet_obj$clustering$clusters)) {
#'     filename <- paste0("heatmap_muscadet_k", k, ".png")
#'     title <- paste(
#'         "Example sample |",
#'         muscadet_obj$clustering$params[["dist_method"]],
#'         muscadet_obj$clustering$params[["hclust_method"]],
#'         "|",
#'         paste0("k=", k) ,
#'         "|",
#'         "Weights of omics:",
#'         paste(muscadet_obj$clustering$params[["weights"]], collapse = ", ")
#'     )
#'     heatmapMuscadet(
#'         muscadet_obj,
#'         k = as.integer(k),
#'         filename = filename,
#'         title = title
#'     )
#' }
#'
heatmapMuscadet <- function(x, filename = NULL, k = NULL, clusters = NULL, title = "",
                            add_bulk_lrr = NULL, show_missing = TRUE, white_scale = c(0.3, 0.7),
                            colors = NULL, quiet = FALSE) {

    # Validate input: x must be a muscadet object
    stopifnot("Input object 'x' must be of class 'muscadet'." = inherits(x, "muscadet"))

    # Validate the muscadet object contains clustering results
    stopifnot(
        "The muscadet object 'x' does not contain clustering data (use clusterMuscadet() to perform clustering of log R ratio data)." =
            !is.null(slot(x, "clustering"))
    )

    # Validate the clustering result for the specified k
    stopifnot(
        "The muscadet object must contain clustering results for the specified k." =
            as.character(k) %in% names(slot(x, "clustering")[["clusters"]])
    )

    # Set to no missing cells if only one omic
    if(length(x@omics) == 1) show_missing <- FALSE

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

    # Default addition of bulk data
    if(is.null(add_bulk_lrr)) {
        add_bulk_lrr <- !is.null(x@bulk.data$log.ratio)
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
    if(!is.null(filename)) stopifnot("The directory doesn't exist" = file.exists(dirname(filename)))

    # Get common and all cells in clustering order
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
        message("Output file: ", ifelse(is.null(filename), "None", filename))
    }

    # Create list of heatmap objects
    list_ht <- lapply(slot(x, "omics"), function(muscomic) {

        pdf(file = NULL)
        ht_opt(message = F)

        # Get chromosome info for features
        coord <- muscomic@coverage$coord.features
        chrom <- factor(coord[coord$keep, "CHROM"], levels = unique(coord[coord$keep, "CHROM"]))

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
            bulk_col <- list(circlize::colorRamp2(c(min(bulk_df$bulk.lrr, na.rm = TRUE),
                                                    median(bulk_df$bulk.lrr, na.rm = TRUE),
                                                    max(bulk_df$bulk.lrr, na.rm = TRUE)),
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
                             row_split = factor(clusters[all_cells], levels=sort(unique(clusters))),
                             row_order = names(clusters),
                             cluster_rows = F, merge_legend = TRUE)


    } else if (show_missing == FALSE) {
        if (is.null(clusters)) {
            # 2. Clustering hclust from the muscadet object clustering with common cells
            hc <- x@clustering$hclust # hclust object to print the dendrogram on the heatmap
            n_cells <- table(dendextend::cutree(hc, k, order_clusters_as_data = FALSE))

            ht_all <- ComplexHeatmap::draw(ht_list, column_title = title, ht_gap = unit(1, "cm"),
                                 cluster_rows = hc, row_split = as.integer(k), row_dend_reorder = FALSE, merge_legend = TRUE)
        } else {
            # 3. Custom cluster assignments vector
            n_cells <- table(clusters)
            ht_all <- ComplexHeatmap::draw(ht_list, column_title = title, ht_gap = unit(1, "cm"),
                                 row_split = factor(clusters[all_cells], levels=sort(unique(clusters))),
                                 row_order = names(clusters), cluster_rows = F, merge_legend = TRUE)
        }
    }
    dev.off()

    # Save complete plot of heatmaps as PNG
    if (!is.null(filename)) {
        png(
            filename = filename,
            width = ht_all@ht_list_param[["width"]],
            height = ht_all@ht_list_param[["height"]],
            units = "mm",
            res = 300
        )
    } else {
        pdf(file = NULL)
    }

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

    # store plot object
    if (is.null(filename)) {
        plot.obj <- list(
            plot = grid.grab(),
            width = ht_all@ht_list_param[["width"]],
            height = ht_all@ht_list_param[["height"]]
        )
    }

    dev.off()

    # return the plot object
    if (is.null(filename)) {
        return(plot.obj)
    }
}





#' Silhouette plot for `muscadet` object
#'
#' Generate a silhouette plot for a specified clustering result within a
#' `muscadet` object.
#'
#' @param x A \code{\link{muscadet}} object containing clustering data (using
#'   [muscadet::clusterMuscadet()]).
#'
#' @param k Integer specifying the number of clusters to plot (`integer`).
#'   It should be within the range of k used for clustering (with
#'   [muscadet::clusterMuscadet()]).
#'
#' @param colors Vector of colors for the cluster annotation (`character`
#'   vector). Default is `NULL`, which uses predefined colors.
#'
#' @param title Character string for the title of the plot (`character`
#'   string). If `NULL`, a default title is generated.
#'
#' @return A ggplot object representing the silhouette plot.
#'
#' @import ggplot2
#' @importFrom stats aggregate
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' library("ggplot2")
#'
#' # Load a muscadet object
#' data(muscadet_obj)
#' plotSil(muscadet_obj, k = 3)
#'
#' # Loop over k
#' for (k in names(muscadet_obj$clustering$clusters)) {
#'     plot <- plotSil(muscadet_obj, k = as.integer(k))
#'     ggsave(paste0("plot_silhouette_", k, ".png"), plot)
#' }
#'
plotSil <- function(x, k, colors = NULL, title = NULL) {

    # Validate input: x must be a muscadet object
    stopifnot("Input object 'x' must be of class 'muscadet'." = inherits(x, "muscadet"))

    # Validate the muscadet object contains clustering results
    stopifnot(
        "The muscadet object 'x' does not contain clustering data (use clusterMuscadet() to perform clustering of log R ratio data)." =
            !is.null(slot(x, "clustering"))
    )

    # Validate the clustering result for the specified k
    stopifnot(
        "The muscadet object must contain clustering results for the specified k." =
            as.character(k) %in% names(slot(x, "clustering")[["clusters"]])
    )

    # Set default color palette for clusters if not provided
    if (is.null(colors)) {
        colors <- c(
            "#FABC2A", "#7FC97F", "#EE6C4D", "#39ADBD", "#BEAED4",
            "#FEE672", "#F76F8E", "#487BEA", "#B67BE6", "#F38D68",
            "#7FD8BE", "#F2AFC0"
        )
    }

    # Generate a default title if none is provided
    if (is.null(title)) {
        title <- paste("Silhouette Plot (k =", as.character(k), ")")
    }

    # Extract silhouette data for the specified k
    sil <- x@clustering[["silhouette"]][["sil.obj"]][[as.character(k)]]
    df <- as.data.frame(sil[, 1:3], stringsAsFactors = TRUE)

    # Order data for plotting
    df <- df[order(df$cluster, -df$sil_width), ]
    df$name <- factor(rownames(df), levels = rownames(df))
    df$cluster <- as.factor(df$cluster)

    # Calculate average silhouette width per cluster
    avg_sil_width <- stats::aggregate(df$sil_width ~ df$cluster, data = df, mean)
    colnames(avg_sil_width) <- c("cluster", "sil_width")

    # Determine positions for cluster annotations
    cluster_counts <- table(df$cluster)
    cluster_positions <- sapply(unique(df$cluster), function(i) {
        # Reverse y-axis positions since it is flipped
        rev_y_position <- sum(cluster_counts) -
            (sum(cluster_counts[as.numeric(levels(df$cluster)[1:i])]) - cluster_counts[i] / 2)
        rev_y_position
    })

    # Create the ggplot object
    p <- ggplot(df, aes(x = .data$sil_width, y = .data$name, fill = .data$cluster)) +
        geom_bar(stat = "identity", width = 1) +
        scale_y_discrete(limits = rev) +
        theme_bw() +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"),
            legend.position = "none"
        ) +
        scale_fill_manual(values = colors) +
        labs(
            x = "Silhouette Width",
            y = "",
            title = title,
            subtitle = paste0(
                "Average Silhouette Width = ", round(mean(df$sil_width), 4), "\n",
                "Annotation: cluster | number of cells | average silhouette width"
            )
        ) +
        ggplot2::xlim(c(min(df$sil_width), max(df$sil_width) + 0.25)) +
        geom_vline(xintercept = mean(df$sil_width), linetype = "dashed", color = "red")

    # Add annotations for each cluster
    for (i in seq_along(cluster_positions)) {
        cluster <- unique(df$cluster)[i]
        avg_width <- round(avg_sil_width$sil_width[avg_sil_width$cluster == cluster], 4)
        count <- cluster_counts[i]
        p <- p + geom_text(
            label = paste(cluster, "|", count, "|", avg_width),
            x = max(df$sil_width) + 0.05,
            y = cluster_positions[i],
            inherit.aes = FALSE,
            hjust = 0,
            vjust = -0.5,
            check_overlap = TRUE
        )
    }

    return(p)
}



#' Plot clustering validation indexes for a `muscadet` object
#'
#' It generates a plot of clustering validation indexes for a clustering result
#' within a `muscadet` object. The index values are computed only using
#' distances between common cells across omics in the `muscadet` object.
#'
#' @param x A \code{\link{muscadet}} object containing clustering data
#'   (generated using [muscadet::clusterMuscadet()]).
#'
#' @param index Character vector specifying one or more validation indexes to
#'   plot among "silhouette", "dunn2", "davisbouldin", "pearsongamma", and "c".
#'   If `NULL`, by default all available indexes are included. If multiple
#'   indexes are selected, the values are normalized for comparability.
#'
#' @param colors Vector of colors for each index in the plot (`character`
#'   vector). Default is `NULL`, which uses predefined colors for the indexes.
#'
#' @param title Character string for the title of the plot (`character`
#'   string). If `NULL`, a default title is generated.
#'
#' @return A ggplot object visualizing the clustering validation indexes across
#'   different numbers of clusters (`k`).
#'
#' @details The function computes several clustering validation indexes,
#'   including:
#'   \itemize{
#'     \item \strong{Silhouette}: Measures how similar an object is to its own
#'     cluster compared to others (see [cluster::silhouette]).
#'     \item \strong{Dunn2}: The ratio of the smallest distance between
#'     observations in different clusters to the largest within-cluster distance
#'     (see [fpc::cluster.stats] `$dunn2`).
#'     \item \strong{Davis-Bouldin}: Measures cluster compactness and separation
#'     (see [clusterSim::index.DB]).
#'     \item \strong{Pearson's Gamma}: Evaluates the goodness of clustering
#'     based on correlation (see [fpc::cluster.stats] `$pearsongamma`).
#'     \item \strong{C Index}: Measures the clustering quality compared to
#'     random data (see [clusterSim::index.C]).
#'   }
#'
#'   If multiple indexes are selected, the values are normalized to fall between
#'   0 and 1. For indexes that are better when minimized ("pearsongamma" and
#'   "c"), their values are reversed for easier comparison.
#'   The k for which the mean of indexes is maximal is highlighted with a dot.
#'
#' @import ggplot2
#' @importFrom stats as.dist
#' @importFrom cluster silhouette
#' @importFrom fpc cluster.stats
#' @importFrom clusterSim index.DB index.C
#' @importFrom tidyr pivot_longer
#'
#' @export
#'
#' @examples
#' # Load muscadet object
#' data(muscadet_obj)
#'
#' # Plot all indexes
#' plotIndexes(muscadet_obj)
#'
#' # Plot specific indexes
#' plotIndexes(muscadet_obj, index = "silhouette")
#'
plotIndexes <- function(x, index = NULL, colors = NULL, title = NULL) {

    # Validate input: x must be a muscadet object
    stopifnot("Input object 'x' must be of class 'muscadet'." = inherits(x, "muscadet"))

    # Validate the muscadet object contains clustering results
    stopifnot(
        "The muscadet object 'x' does not contain clustering data (use clusterMuscadet() to perform clustering of log R ratio data)." =
            !is.null(slot(x, "clustering"))
    )

    # Define default indexes if none are provided
    if (is.null(index)) {
        index <- c("silhouette", "dunn2", "davisbouldin", "pearsongamma", "c")
    }
    # Check for indexes correct names
    index <- match.arg(index, c("silhouette", "dunn2", "davisbouldin", "pearsongamma", "c"), several.ok = TRUE)

    # Set default colors if not provided
    if (is.null(colors)) {
        colors <- c(
            "silhouette" = "brown",
            "dunn2" = "coral2",
            "davisbouldin" = "tan2",
            "pearsongamma" = "turquoise4",
            "c" = "skyblue4"
        )
    }

    # Generate a default title if none is provided
    if (is.null(title)) {
        title <- "Clustering Validation Indexes Across k"
    }

    # Extract the clustering data
    k_clusters <- as.integer(names(x@clustering$clusters))  # Number of clusters (k)
    dist <- stats::as.dist(x@clustering$dist)  # Distance matrix to dist object

    # Initialize a data frame to store index values
    df_indexes <- data.frame(k = k_clusters)

    # Compute selected indexes for each k
    for (k in k_clusters) {

        # Extract cluster assignments for the current k
        clusters <- x@clustering$clusters[[as.character(k)]]
        clusters <- clusters[names(dist)] # Restrict to common cells in dist

        # Compute each index if selected
        if ("silhouette" %in% index) {
            df_indexes[df_indexes$k == k, "silhouette"] <-
                summary(cluster::silhouette(as.integer(clusters), dist))[["avg.width"]]
        }
        if ("dunn2" %in% index) {
            df_indexes[df_indexes$k == k, "dunn2"] <-
                fpc::cluster.stats(dist, as.integer(clusters))$dunn2
        }
        if ("davisbouldin" %in% index) {
            df_indexes[df_indexes$k == k, "davisbouldin"] <-
                clusterSim::index.DB(dist, as.integer(clusters))$DB
        }
        if ("pearsongamma" %in% index) {
            df_indexes[df_indexes$k == k, "pearsongamma"] <-
                fpc::cluster.stats(dist, as.integer(clusters))$pearsongamma
        }
        if ("c" %in% index) {
            df_indexes[df_indexes$k == k, "c"] <-
                clusterSim::index.C(dist, as.integer(clusters))
        }
    }

    # If multiple indexes selected:
    if( length(index) > 1) {
        # Normalize indexes to values between 0 and 1
        df_indexes[, 2:ncol(df_indexes)] <- apply(df_indexes[, 2:ncol(df_indexes)], 2, function(x) {
            (x - min(x)) / diff(range(x))
        })
        # Reverse indexes that should be minimized, for comparability
        if ("pearsongamma" %in% index) {
            df_indexes <- mutate(df_indexes, pearsongamma = 1 - .data$pearsongamma)
        }
        if ("c" %in% index) {
            df_indexes <- mutate(df_indexes, c = 1 - .data$c)
        }
    }

    # Find optimal k based on selected indexes
    if (ncol(df_indexes) > 2) {
        # Find the k with the maximum mean index value
        opt_k <- df_indexes[which(rowMeans(df_indexes[, 2:ncol(df_indexes)]) == max(rowMeans(df_indexes[, 2:ncol(df_indexes)]))), "k"]
    } else {
        if (colnames(df_indexes)[2] %in% c("pearsongamma", "c")) {
            # Find k with the minimal value for indexes to minimize
            opt_k <- df_indexes[which(df_indexes[, 2] == min(df_indexes[, 2])), "k"]
        } else {
            # Find k with the maximal value for indexes to maximize
            opt_k <- df_indexes[which(df_indexes[, 2] == max(df_indexes[, 2])), "k"]
        }
    }

    # Transform the data frame for plotting
    df_plot <- tidyr::pivot_longer(df_indexes, -k, names_to = "Index", values_to = "Value")
    df_plot$Index <- factor(df_plot$Index, levels = index) # Maintain order of indexes

    # Generate the plot
    plot <- ggplot(df_plot, aes(
        x = .data$k,
        y = .data$Value,
        color = .data$Index,
        group = .data$Index
    )) +
        geom_line(linewidth = 1) +
        geom_point(data = df_plot[which(df_plot$k == opt_k),]) +
        labs(x = "Number of Clusters (k)",
             y = "Normalized Index Value",
             title = title) +
        scale_x_continuous(breaks = unique(df_plot$k)) +
        scale_color_manual(
            values = colors,
            labels = c(
                "silhouette" = "Silhouette",
                "dunn2" = "Dunn2",
                "pearsongamma" = "Pearson's Gamma",
                "davisbouldin" = "Davis-Bouldin",
                "c" = "C Index"
            )
        ) +
        theme_classic() +
        theme(legend.position = "right")

    return(plot)
}
