% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dotplot_gsea.R
\name{dotplot_gsea}
\alias{dotplot_gsea}
\title{Dotplot for fgsea_all output with normalized enrichment scores and dots sized
by adjusted p-values}
\usage{
dotplot_gsea(fgsea_all_res, top = 7, pad = 0.1, wrap = 50, ...)
}
\arguments{
\item{fgsea_all_res}{a tibble from fgsea_all}

\item{top}{Number of top pathways by NES to display, default 7}

\item{pad}{padding for x-axis, default 0.1}

\item{wrap}{wrap pathway strings, default 50 characters}

\item{\dots}{additional options passed to scale_x_continuous}
}
\value{
a ggplot
}
\description{
Dotplot for fgsea_all output with normalized enrichment scores and dots sized
by adjusted p-values
}
\examples{
\dontrun{
h1 <- fgsea_all(res, msig_hallmark)
dotplot_gsea(h1[[1]])
# pass labels to scale_x_continuous if too many digits
dotplot_gsea(h1[[1]], labels = function(x) format(round(x, digits=2), nsmall = 2))
}
}
\author{
Chris Stubben
}
