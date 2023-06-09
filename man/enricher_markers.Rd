% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enricher_markers.R
\name{enricher_markers}
\alias{enricher_markers}
\title{Run enricher from clusterProfiler on findAllMarkers}
\usage{
enricher_markers(res, gsets, maxGSSize = 2000, padj = 0.05, quiet = TRUE, ...)
}
\arguments{
\item{res}{A list of DESeq results}

\item{gsets}{A list of gene sets (or two column dataframe with pathway and gene)}

\item{maxGSSize}{Maximum set size, default 2000}

\item{padj}{Adjusted p-value cutoff for significant genes, default 0.05}

\item{quiet}{suppress messages except total significant}

\item{\dots}{Additional options like minGSSize and pvalueCutoff passed to \code{enricher}}
}
\value{
A list of tibbles
}
\description{
Run enricher from clusterProfiler on findAllMarkers
}
\examples{
\dontrun{
library(hciRdata)
enricher_all(res, msig_pathways$KEGG)
}
}
\author{
Chris Stubben
}
