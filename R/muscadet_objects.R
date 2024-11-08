# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions ------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The muscomic class
#'
#' A `muscomic` object encapsulates coverage and allelic data from one
#' single-cell omic dataset used as primary input for muscadet analysis.
#'
#' @slot type Type of single-cell omic. Either ATAC or RNA currently supported
#'   (`character`).
#' @slot label.omic Label to display for the single-cell omic (`character`).
#' @slot coverage Coverage data on features (`list`).
#' @slot allelic Allelic data at single nucleotide polymorphisms (SNPs) positions (`list`).
#'
#' @aliases muscomic
#'
#' @importFrom methods setClass
#' @exportClass muscomic

methods::setClass(
  "muscomic",
  slots = c(
    type = "character",
    label.omic = "character",
    coverage = "list",
    allelic = "list"
  )
)


#' The muscadet class
#'
#' A `muscadet` object encapsulates data from different single-cell omics
#' as \code{\link{muscomic}} objects as well as downstream analysis results after
#' clustering and CNA calling.
#'
#' @slot omics List of \code{\link{muscomic}} objects, one per single-cell omic (`list`).
#' @slot bulk.data List of objects containing data from paired bulk sequencing
#'   (`list`).
#' @slot clustering List of objects containing data associated with the
#'   clustering of cells based on coverage log R ratio values (`list`).
#' @slot cnacalling List of objects containing data associated with the calling
#'   of copy number alterations (CNAs) (`list`).
#' @slot genome Reference genome name among: "hg38", "hg19" and "mm10" (`character`).
#'
#' @aliases muscadet
#'
#' @importFrom methods setClass
#' @exportClass muscadet

methods::setClass(
  "muscadet",
  slots = c(
    omics = "list",
    bulk.data = "list",
    clustering = "list",
    cnacalling = "list",
    genome = "character"
  )
)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create objects functions -----------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create a muscomic object
#'
#' Create a \code{\link{muscomic}} object.
#'
#' @param type.omic Type of single cell omic, either "ATAC" or "RNA"
#'   (`character` string).
#' @param mat_counts Matrix of raw counts *cells x features* (`matrix` or
#'   `dgCMatrix`).
#' @param allele_counts Data frame of allele counts at single nucleotide
#'   polymorphisms (SNPs) positions per cell (`data.frame`).
#' @param features Data frame of features (peaks, genes...) coordinates on
#'   genome (`data.frame`).
#' @param label.omic Label for the single cell omic (`character` string). By
#'   default "scATAC-seq" is used for "ATAC" type and "scRNA-seq" for "RNA"
#'   type.
#' @param label.features Label for features (`character` string). By default
#'   "peaks" is used for "ATAC" type and "genes" for "RNA" type.
#'
#' @return
#' A \code{\link{muscomic}} object.
#'
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom Matrix Matrix
#' @importFrom stringr str_remove
#' @importFrom gtools mixedsort
#'
#' @export
#'
#' @seealso
#' * [CreateMuscadetObject()]
#' * \code{\link{muscomic-class}}, \code{\link{muscadet-class}}
#'
#' @examples
#' atac <- CreateMuscomicObject(
#'   type = "ATAC",
#'   mat_counts = mat_counts_atac_tumor,
#'   allele_counts = allele_counts_atac_tumor,
#'   features = peaks
#' )
#'
#' rna <- CreateMuscomicObject(
#'   type = "RNA",
#'   mat_counts = mat_counts_rna_tumor,
#'   allele_counts = allele_counts_rna_tumor,
#'   features = genes
#' )
#'
#' atac_ref <- CreateMuscomicObject(
#'   type = "ATAC",
#'   mat_counts = mat_counts_atac_ref,
#'   allele_counts = allele_counts_atac_ref,
#'   features = peaks
#' )
#'
#' rna_ref <- CreateMuscomicObject(
#'   type = "RNA",
#'   mat_counts = mat_counts_rna_ref,
#'   allele_counts = allele_counts_rna_ref,
#'   features = genes
#' )
CreateMuscomicObject <- function(type.omic = c("ATAC", "RNA"),
                                 mat_counts,
                                 allele_counts,
                                 features,
                                 label.omic = NULL,
                                 label.features = NULL) {
  # check for type of omic
  type.omic <- match.arg(arg = type.omic)

  # default labels
  if (type.omic == "ATAC") {
    if (is.null(label.omic)) {
      label.omic <- "scATAC-seq"
    }
    if (is.null(label.features)) {
      label.features <- "peaks"
    }
  }
  if (type.omic == "RNA") {
    if (is.null(label.omic)) {
      label.omic <- "scRNA-seq"
    }
    if (is.null(label.features)) {
      label.features <- "genes"
    }
  }

  # check count matrix
  stopifnot(
    "`mat_counts` must be either a matrix or a dgCMatrix." =
      any(class(mat_counts) %in% c("matrix", "dgCMatrix"))
  )
  # convert to dgCMatrix class objects if matrix
  if (any(class(mat_counts) == "matrix")) {
    mat_counts <- Matrix::Matrix(mat_counts)
  }

  # check coordinates of features
  stopifnot("`features` must be a data frame." = class(features) == "data.frame")
  features$CHROM <- stringr::str_remove(features$CHROM, "chr") # remove "chr" if necessary
  features$CHROM <- ordered(features$CHROM,
    levels = gtools::mixedsort(unique(features$CHROM))
  ) # ordered chromosomes
  features <- dplyr::arrange(features, CHROM, start) # sort by chromosome and position

  # check number of features is identical in coordinates data frame and in number of rows of matrix
  stopifnot(
    "Feature coordinates data frame and count matrix must have the same number of rows." =
      nrow(features) == nrow(mat_counts)
  )

  # check matching cells between coverage (count matrix) and allelic data (allele count df)
  stopifnot(
    "Cells in the allele counts data frame must be the same as the ones in matrix counts column names." =
      unique(allele_counts$cell) %in% colnames(mat_counts)
  )

  # generate table of raw counts by converting matrix to summary data frame
  coord.df <- dplyr::mutate(features, index = as.numeric(1:nrow(features))) # unique index value

  table_counts <- as.data.frame(Matrix::summary(mat_counts)) %>%
    dplyr::mutate(cell = colnames(mat_counts)[j], omic = type.omic) %>%
    dplyr::rename(index = i, DP = x) %>%
    dplyr::left_join(coord.df, by = "index") %>%
    dplyr::select(c(omic, cell, id, CHROM, start, end, DP))

  # coverage data
  coverage <- list(
    mat.counts = mat_counts,
    table.counts = table_counts,
    coord.features = features,
    label.features = label.features
  )
  # allelic data
  allelic <- list(table.counts = data.frame(omic = type.omic, allele_counts))

  # create object
  obj <- new(
    Class = "muscomic",
    type = type.omic,
    label.omic = label.omic,
    coverage = coverage,
    allelic = allelic
  )
  return(obj)
}


#' Create a muscadet object
#'
#' Create a \code{\link{muscadet}} object.
#'
#' @param omics A list of \code{\link{muscomic}} objects (`list`). The names of
#'   the list will set the names of omics in the final object, if the list is
#'   unamed, the type is taken instead.
#' @param bulk_lrr A data frame containing log R ratio per genomic segments from
#'   bulk sequencing data (`data.frame`).
#' @param bulk.label Label for bulk data (`character` string).
#' @param genome Reference genome name among: "hg38", "hg19" and "mm10"
#'   (`character` string). "hg38" by default.
#'
#' @return
#' A \code{\link{muscadet}} object.
#'
#' @export
#'
#' @note A \code{\link{muscadet}} object can contain several
#'   \code{\link{muscomic}} objects of the same type (`ATAC` or `RNA` for slot
#'   `type`) but they can't have identical labels (slot `label.omic`).
#'
#' @seealso
#' * [CreateMuscomicObject()]
#' * \code{\link{muscomic-class}}, \code{\link{muscadet-class}}
#'
#' @examples
#' muscadet <- CreateMuscadetObject(
#'   omics = list(atac, rna),
#'   bulk_lrr = bulk_lrr,
#'   bulk.label = "WGS",
#'   genome = "hg38"
#' )
#' muscadet_ref <- CreateMuscadetObject(
#'   omics = list(atac_ref, rna_ref),
#'   genome = "hg38"
#' )
#'
CreateMuscadetObject <- function(omics,
                                 bulk_lrr = NULL,
                                 bulk.label = NULL,
                                 genome = "hg38") {
  # check for omic objects
  for (i in seq_along(omics)) {
    stopifnot(
      "`omics` list elements must be of class `muscomic` (use CreateMuscomicObject())." =
        as.character(class(omics[[i]])) == "muscomic"
    )
  }


  # set names of omics based on type if unamed
  if (is.null(names(omics))) {
    omics <- setNames(omics, unlist(lapply(omics, function(o) {
      o@type
    })))

    omics <- setNames(omics, ave(
      names(omics),
      names(omics),
      FUN = function(i) {
        if (length(i) > 1) {
          paste0(i, "_", seq_along(i))
        } else {
          i
        }
      }
    ))
  }

  # labels of omics can't be identical
  stopifnot(
    "Identical omic labels found in the `omics` (muscomic object list) provided. You can check the labels with `slot(muscomic_oj, 'label.omic')`." =
      !any(duplicated(unlist(
        lapply(omics, function(o) {
          o@label.omic
        })
      )))
  )

  # check for genome
  stopifnot(
    "Genome must be either 'hg38', 'hg19' or 'mm10'." =
      genome %in% c("hg38", "hg19", "mm10")
  )

  # add bulk log R ratio data
  bulk.data <- list(log.ratio = bulk_lrr, label = bulk.label)

  # create object
  obj <- new(
    Class = "muscadet",
    omics = omics,
    bulk.data = bulk.data,
    genome = genome
  )
  return(obj)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods -------------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# show -------------------------------------------------------------------------

#' muscomic object overview
#'
#' Overview of a \code{\link{muscomic}} object.
#'
#' @param object A \code{\link{muscomic}} object.
#'
#' @return Prints summary to \code{\link[base]{stdout}} and invisibly returns
#'   \code{NULL}
#'
#' @keywords internal
#'
#' @seealso \code{\link{muscomic-class}}, [CreateMuscomicObject()]
#'
#' @examples
#' atac <- CreateMuscomicObject(
#'   type = "ATAC",
#'   mat_counts = mat_counts_atac_tumor,
#'   allele_counts = allele_counts_atac_tumor,
#'   features = peaks
#' )
#' atac
setMethod(
  f = "show",
  signature = signature(object = "muscomic"),
  definition = function(object) {
    cat(
      "A muscomic object of type",
      object@type,
      "labelled", object@label.omic, "containing:",
      "\n",
      ncol(object@coverage$mat.counts),
      "cells",
      "\n",
      nrow(object@coverage$mat.counts),
      object@coverage$label.features,
      "\n",
      length(unique(atac@allelic$table.counts$id)),
      "SNPs",
      "\n"
    )
  }
)

#' muscadet object overview
#'
#' Overview of a \code{\link{muscadet}} object.
#'
#' @param object A \code{\link{muscadet}} object.
#'
#' @return Prints summary to \code{\link[base]{stdout}} and invisibly returns
#'   \code{NULL}
#'
#' @keywords internal
#'
#' @seealso \code{\link{muscadet-class}}, [CreateMuscadetObject()]
#'
#' @examples
#' atac <- CreateMuscomicObject(
#'   type = "ATAC",
#'   mat_counts = mat_counts_atac_tumor,
#'   allele_counts = allele_counts_atac_tumor,
#'   features = peaks
#' )
#' rna <- CreateMuscomicObject(
#'   type = "RNA",
#'   mat_counts = mat_counts_rna_tumor,
#'   allele_counts = allele_counts_rna_tumor,
#'   features = genes
#' )
#'
#' muscadet <- CreateMuscadetObject(
#'   omics = list(atac, rna),
#'   bulk_lrr = bulk_lrr,
#'   bulk.label = "WGS",
#'   genome = "hg38"
#' )
setMethod(
  f = "show",
  signature = signature(object = "muscadet"),
  definition = function(object) {
    cat(
      "A muscadet object",
      "\n",
      "omics:", paste(names(object@omics), collapse = ", "), "\n",
      " ", "types:", paste(lapply(object@omics, function(i) {
        i@type
      }), collapse = ", "), "\n",
      " ", "labels:", paste(lapply(object@omics, function(i) {
        i@label.omic
      }), collapse = ", "), "\n",
      " ", "cells:", paste(lapply(object@omics, function(i) {
        ncol(i@coverage$mat.counts)
      }), collapse = ", "), "\n",
      " ", "features:", paste(lapply(object@omics, function(i) {
        nrow(i@coverage$mat.counts)
      }), collapse = ", "), "\n",
      " ", "feature labels:", paste(lapply(object@omics, function(i) {
        i@coverage$label.features
      }), collapse = ", "), "\n",
      " ", "SNPs:", paste(lapply(object@omics, function(i) {
        length(unique(i@allelic$table.counts$id))
      }), collapse = ", "), "\n",
      "data from paired bulk sequencing:", ifelse(is.null(object@bulk.data$label), "None", object@bulk.data$label), "\n",
      "clustering:", ifelse(length(object@clustering) == 0, "None", object@clustering), "\n",
      "CNA calling:", ifelse(length(object@cnacalling) == 0, "None", object@cnacalling), "\n",
      "genome:", object@genome
    )
  }
)



# New methods ------------------------------------------------------------------

#' \code{\link{muscomic}} and \code{\link{muscadet}} methods
#'
#' Methods for \code{\link{muscomic}} and \code{\link{muscadet}} objects, to
#' facilitate access to slots within the objects.
#'
#' @name muscadet-methods
#' @rdname muscadet-methods
#' @aliases muscomic-methods
#'
#' @param x A \code{\link{muscomic}} or \code{\link{muscadet}} object.
#'
#' @seealso \code{\link{muscomic-class}} \code{\link{muscadet-class}}
#'
NULL

#' @rdname muscadet-methods
#'
#' @return
#' `Cells`:
#' - `muscomic`: A vector of cell names.
#' - `muscadet`: A list of vectors of cell names, one list element per omic.
#'
Cells <- function(x, ...) {
  UseMethod(generic = "Cells", object = x)
}

#' @rdname muscadet-methods
#'
#' @return
#' `Features`:
#' - `muscomic`: A vector of feature names.
#' - `muscadet`: A list of vectors of feature names, one list element per omic.
#'
Features <- function(x, ...) {
  UseMethod(generic = "Features", object = x)
}

#' @rdname muscadet-methods
#'
#' @return
#' `coordFeatures`:
#' - `muscomic`: A data frame of feature coordinates.
#' - `muscadet`: A list of data frames of feature coordinates, one list element per omic.
#'
coordFeatures <- function(x, ...) {
  UseMethod(generic = "coordFeatures", object = x)
}

#' @rdname muscadet-methods
#'
#' @return
#' `matCounts`:
#' - `muscomic`: A \code{\link{dgCMatrix}} *features x cells*.
#' - `muscadet`: A list of \code{\link{dgCMatrix}} *features x cells*, one list element per omic.
#'
matCounts <- function(x, ...) {
  UseMethod(generic = "matCounts", object = x)
}

#' @rdname muscadet-methods
#'
setMethod(
  f = "Cells",
  signature = signature(x = "muscomic"),
  definition = function(x) {
    return(colnames(slot(x, "coverage")[["mat.counts"]]))
  }
)

#' @rdname muscadet-methods
#'
setMethod(
  f = "Cells",
  signature = signature(x = "muscadet"),
  definition = function(x) {
    lapply(slot(x, "omics"), function(omic) {
      return(colnames(slot(omic, "coverage")[["mat.counts"]]))
    })
  }
)

#' @rdname muscadet-methods
#'
setMethod(
  f = "Features",
  signature = signature(x = "muscomic"),
  definition = function(x) {
    return(rownames(slot(x, "coverage")[["mat.counts"]]))
  }
)

#' @rdname muscadet-methods
#'
setMethod(
  f = "Features",
  signature = signature(x = "muscadet"),
  definition = function(x) {
    lapply(slot(x, "omics"), function(omic) {
      return(rownames(slot(omic, "coverage")[["mat.counts"]]))
    })
  }
)

#' @rdname muscadet-methods
#'
setMethod(
  f = "coordFeatures",
  signature = signature(x = "muscomic"),
  definition = function(x) {
    return(slot(x, "coverage")[["coord.features"]])
  }
)

#' @rdname muscadet-methods
#'
setMethod(
  f = "coordFeatures",
  signature = signature(x = "muscadet"),
  definition = function(x) {
    lapply(slot(x, "omics"), function(omic) {
      return(slot(omic, "coverage")[["coord.features"]])
    })
  }
)

#' @rdname muscadet-methods
#'
setMethod(
  f = "matCounts",
  signature = signature(x = "muscomic"),
  definition = function(x) {
    return(slot(x, "coverage")[["mat.counts"]])
  }
)

#' @rdname muscadet-methods
#'
setMethod(
  f = "matCounts",
  signature = signature(x = "muscadet"),
  definition = function(x) {
    lapply(slot(x, "omics"), function(omic) {
      return(slot(omic, "coverage")[["mat.counts"]])
    })
  }
)
