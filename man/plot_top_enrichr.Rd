% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_top_enrichr.R
\name{plot_top_enrichr}
\alias{plot_top_enrichr}
\title{Cluster -log10 pvalues from top gene sets in each contrast}
\usage{
plot_top_enrichr(
  e1,
  db = 1,
  top_n = 5,
  trim = 50,
  min_p = 100,
  complexHeatmap = TRUE,
  ...
)
}
\arguments{
\item{e1}{results from enrichr_all}

\item{db}{List index or name, default 1}

\item{top_n}{Number of top gene sets to plot}

\item{trim}{Trim long names, default 50}

\item{min_p}{Minimum p-value to plot (to adjust scale)}

\item{complexHeatmap}{use pheatmap from ComplexHeatmap}

\item{\dots}{additional options passed to pheatmap}
}
\value{
pheatmp
}
\description{
Cluster -log10 pvalues from top gene sets in each contrast
}
\examples{
\dontrun{
library(enrichR)
setEnrichrSite("Enrichr")
celldbs <- c("KEGG_2021_Human", "WikiPathways_2019_Human")
e1 <- enrichr_all(res, celldbs)
plot_top_enrichr(e1, 1, top=3)
}
}
\author{
Chris Stubben
}
