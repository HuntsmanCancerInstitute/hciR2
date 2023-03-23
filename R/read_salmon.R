#' Read Salmon output files
#'
#' Reads Salmon counts or stats files and optionally reshape into wide format.
#'
#' @param path the path to Salmon output files, the default corresponds to the working directory.
#' @param pattern regular expression for file name matching, default .sf
#' @param TPM read TPM counts, default is NumReads
#' @param reshape reshape into wide format with samples in rows (a count matrix).
#'
#' @return A tibble in long or wide format if reshape=TRUE
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  # raw counts
#'   counts <- read_salmon( "Alignments")
#'  # TPMs
#'  tpm <- read_salmon("Alignments", TPM=TRUE)
#' }
#' @export

read_salmon <- function(path = ".", pattern, TPM = FALSE, reshape = TRUE){
     if(missing(pattern))  pattern <- "\\.sf$"
     fc <- read_sample_files(path, pattern)
    if(reshape){
       if(TPM){
         fc  <- dplyr::select(fc, sample, Name, TPM) %>% tidyr::spread(sample, TPM)
       }else{
         fc  <- dplyr::select(fc, sample, Name, NumReads) %>% tidyr::spread(sample, NumReads)
      }
   fc  <-  fc[, c(1, order_samples(colnames(fc)[-1])+1) ]
    }
  fc
}
