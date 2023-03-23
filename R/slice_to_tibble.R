#' Slice a coverage vector and convert to tibble
#'
#' Converts the IRanges output from slice to a tibble with peak widths and
#' heights
#'
#' @param cov Coverage RleList
#' @param lower lower bound of slice
#'
#' @return A tibble with peak width and height
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'   rbga <- readGAlignments("14415X1.bam")
#'   bamCov <- coverage(rbga)
#'   slice_to_tibble(bamCov, lower = 100)
#' }
#' @export

slice_to_tibble <- function( cov, lower= 100 ){
  ## IRanges with start, end, width
  slices <- unlist( IRanges::slice(cov, lower = lower, rangesOnly=TRUE) )
  ## data.frame with start, end, width, names (chr)
  x <- data.frame(slices)
  x <- x[, c(4, 1:3)]
  names(x)[1] <- "seqnames"
  # ADD height
  x$height <- unlist( viewMaxs( IRanges::slice(cov, lower = lower) ))
  slices <- x %>% arrange(desc(height)) %>% as_tibble()
  message("Found ", nrow(slices), " peaks")
  slices
}
