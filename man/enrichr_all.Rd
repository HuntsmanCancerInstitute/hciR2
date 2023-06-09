% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichr_all.R
\name{enrichr_all}
\alias{enrichr_all}
\title{#' Run Enrichr on all DESeq2 result tables}
\usage{
enrichr_all(
  res,
  celldbs,
  deseq.padj = 0.05,
  fisher.padj = 0.05,
  logFC,
  protein_coding = TRUE
)
}
\arguments{
\item{res}{A list of DESeq results}

\item{celldbs}{A vector of gene set names from listEnrichrDbs}

\item{deseq.padj}{cutoff for significant genes to test, default 0.05}

\item{fisher.padj}{cutoff for testing overlap, default 0.05}

\item{logFC}{log2FoldChange cutoff for significant genes, default none}

\item{protein_coding}{Test protein coding genes only, default TRUE}
}
\value{
tibble
}
\description{
#' Run Enrichr on all DESeq2 result tables
}
\note{
if Enrichr is not responding, check SSL certificate
 httr::set_config(httr::config(ssl_verifypeer = 0L))
}
\examples{
\dontrun{
library(enrichR)
setEnrichrSite("Enrichr")
celldbs <- c("KEGG_2021_Human", "WikiPathways_2019_Human")
e1 <- enrichr_all(res, celldbs)
}
}
\author{
Chris Stubben
}
