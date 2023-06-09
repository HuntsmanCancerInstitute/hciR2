% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_salmon.R
\name{read_salmon}
\alias{read_salmon}
\title{Read Salmon output files}
\usage{
read_salmon(path = ".", pattern, TPM = FALSE, reshape = TRUE)
}
\arguments{
\item{path}{the path to Salmon output files, the default corresponds to the working directory.}

\item{pattern}{regular expression for file name matching, default .sf}

\item{TPM}{read TPM counts, default is NumReads}

\item{reshape}{reshape into wide format with samples in rows (a count matrix).}
}
\value{
A tibble in long or wide format if reshape=TRUE
}
\description{
Reads Salmon counts or stats files and optionally reshape into wide format.
}
\examples{
\dontrun{
 # raw counts
  counts <- read_salmon( "Alignments")
 # TPMs
 tpm <- read_salmon("Alignments", TPM=TRUE)
}
}
\author{
Chris Stubben
}
