% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_deseq_results.R
\name{write_deseq_results}
\alias{write_deseq_results}
\title{Write DESeq results}
\usage{
write_deseq_results(result_all, sets, samples, file = "DESeq.xlsx")
}
\arguments{
\item{result_all}{a list from \code{results_all}}

\item{sets}{table with set intersections, optional}

\item{file}{file name}
}
\value{
A tibble
}
\description{
Write DESeq result files to an Excel file.
}
\examples{
\dontrun{
 write_deseq_results(res)
}
}
\author{
Chris Stubben
}
