% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plotProfile}
\alias{plotProfile}
\title{Plot CNA profiles from muscadet object}
\usage{
plotProfile(
  x,
  data,
  title = NULL,
  allelic.type = "lor",
  point.size = 2,
  chrom.colors = c("slategrey", "skyblue"),
  lor.colors = c("peachpuff2", "paleturquoise3"),
  cn.colors = c("black", "brown2"),
  cf.colors = c("white", "steelblue", "bisque2"),
  dipLogR.color = c("magenta4"),
  seg.color = c("brown2")
)
}
\arguments{
\item{x}{A \code{\link{muscadet}} object containing CNA calling data to be
visualized (generated using \code{\link[=cnaCalling]{cnaCalling()}}).}

\item{data}{Either a cluster identifier to plot data of a cluster or
"allcells" to plot data on all cells.}

\item{title}{An optional title for the plot. Default is \code{NULL}.}

\item{allelic.type}{A character string indicating the allelic metric to plot:
"lor" for log odds ratio or "vaf" for variant allele frequency. Default is
"lor".}

\item{point.size}{Numeric value specifying the size of points in the plot in
pixel (with pch = "."). Default is \code{2}.}

\item{chrom.colors}{A character vector of length 2 defining alternating
chromosome colors. Default is \code{c("slategrey", "skyblue")}.}

\item{lor.colors}{A character vector of length 2 for log odds ratio point
colors depending of variant allele frequency in all cells. Use "none" to
use the alternating chromosome colors (defined by \code{chrom.colors}). Default
is \code{c("peachpuff2", "paleturquoise3")}.}

\item{cn.colors}{A character vector of length 2 for total copy number and
minor allele copy number segment colors. Default is \code{c("black", "brown2")}.}

\item{cf.colors}{A character vector of length 3 for cellular fraction
gradient (of 10 values): start color of the gradient, end color of the
gradient, and color for normal diploid (depending on the ploidy). Default
is \code{c("white", "steelblue", "bisque2")}.}

\item{dipLogR.color}{A character string for the diploid log R ratio line
color. Default is "magenta4".}

\item{seg.color}{A character string for the color of segment medians. Default
is "brown2".}
}
\value{
A multi-panel plot of CNA profiles is produced.
}
\description{
This function generates a multi-panel plot of copy number alteration (CNA)
profiles from a \code{\link{muscadet}} object, including: log R ratios
values, log odds ratio (or variant allele frequency), copy numbers and cell
fractions.
}
\examples{
# Load example muscadet object
data(muscadet_obj)

# Plot profile for all cells
pdf("CNAprofile_allcells.pdf", width = 15, height = 7.5) # Save as PDF
plotProfile(muscadet_obj, data = "allcells", title = "Example data - all cells")
dev.off()

}
