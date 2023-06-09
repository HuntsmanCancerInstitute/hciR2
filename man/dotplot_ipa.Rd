% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dotplot_ipa.R
\name{dotplot_ipa}
\alias{dotplot_ipa}
\title{Dotplot for IPA output with normalized enrichment scores and dots sized
by adjusted p-values}
\usage{
dotplot_ipa(
  ipa,
  top = 7,
  pvalue = 0.05,
  wrap = 30,
  pad = 0.1,
  first = "inhibited",
  ...
)
}
\arguments{
\item{ipa}{a table from read_IPA}

\item{top}{Number of top pathways, regulators or functions to display, default 7}

\item{pvalue}{P-value cutoff, default 0.05}

\item{wrap}{wrap pathway strings, default 40 characters}

\item{pad}{padding for x-axis, default 0.1}

\item{first}{plot inhibited values first (down diagonal)}

\item{\dots}{additional options passed to scale_x_continuous}
}
\value{
a ggplot
}
\description{
Dotplot for IPA output with normalized enrichment scores and dots sized
by adjusted p-values
}
\examples{
\dontrun{
ipa1 <- read_ipa("IPA_export_all.txt")
dotplot_ipa(ipa1)
}
}
\author{
Chris Stubben
}
