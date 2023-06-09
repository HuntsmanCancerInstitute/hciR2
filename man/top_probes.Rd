% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/top_probes.R
\name{top_probes}
\alias{top_probes}
\title{Get a count matrix of top probes}
\usage{
top_probes(
  res,
  eset,
  padj = 0.05,
  n = 40,
  row_names = "Gene.Symbol",
  by = "probe"
)
}
\arguments{
\item{res}{results from top_tibble}

\item{eset}{expression set with normalized counts}

\item{padj}{adjusted p-value cutoff, default 0.05}

\item{n}{number of top probes, default 40}

\item{row_names}{row names, Gene.Symbol default}

\item{by}{join results to count rownames by this column, default probe}
}
\value{
A tibble with colData attribute
}
\description{
Join DESeq results to rlog or other counts, mainly for plotting gene heatmaps.
}
\examples{
\dontrun{
 x <- top_counts(res[[1]], eset)
 pheatmap(as_matrix(x))
 plot_genes(x, "trt")
}
}
\author{
Chris Stubben
}
