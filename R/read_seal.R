#' Read SEAL count files and merge into a single matrix
#'
#' @param path the path to SEAL output files, the default
#'        corresponds to the working directory.
#' @param pattern regular expression for file name matching, default .txt
#' @param reshape reshape into wide format with samples in rows (count matrix)
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'    read_seal()
#' }
#' @export

read_seal <- function(path=".", pattern, reshape=TRUE){
     if(missing(pattern))  pattern <- "txt$"
     c1 <- read_sample_files( path, pattern, skip=3)
     names(c1)[1:2] <- c("id", "name")
     c1$sgRNA <- gsub(".*_", "", c1$name)
     c1$Gene <- gsub("_.*", "", c1$name)

     x <- dplyr::select(c1, id, sgRNA, Gene, Reads)
     if(reshape){
         x <- tidyr::spread(x,  id, Reads, fill=0)
         if(ncol(x) > 12) x <- x[, c(1, 2, order_samples(colnames(x)[-c(1:2)])+1) ]
     }
   x
}
