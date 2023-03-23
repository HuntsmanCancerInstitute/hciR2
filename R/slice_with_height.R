#' Slice a stranded genome alignment and add peak height
#'
#' Converts slice output to a GRanges object and adds height using viewMaxs
#'
#' @param aln a GAlignments object
#' @param lower lower bound of slice
#' @param min_width remove peaks with 10 or fewer bases
#'
#' @return A GRanges
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'   rbga <- readGAlignments("14415X1.bam")
#'   slice_with_height(rbga, lower = 100)
#' }
#' @export


slice_with_height <- function( aln, lower= 100, min_width=10){
  ## coverage by strand (add option for unstranded)
  c1p <- IRanges::coverage(subset(aln, strand=="+"))
  c1n <- IRanges::coverage(subset(aln, strand=="-"))
  ## IRanges with start, end, width
  s1 <- unlist( IRanges::slice(c1p, lower = lower, rangesOnly=TRUE) )
  s2 <- unlist( IRanges::slice(c1n, lower = lower, rangesOnly=TRUE) )
  if(length(s1)==0){
    gr1 <- GRanges()
  }else{
    h1 <- unlist( IRanges::viewMaxs( IRanges::slice(c1p, lower = lower) ))
    gr1 <- GRanges(seqnames = names(s1), s1, strand="+",
     length=width(s1), height = h1)
  }
  if(length(s2)==0){
    gr2 <- GRanges()
  }else{
    h2 <- unlist( IRanges::viewMaxs( IRanges::slice(c1n, lower = lower) ))
    gr2 <- GRanges(seqnames = names(s2), s2, strand="-",
     length=width(s2), height = h2)
  }
  ## avoid warning about sequence levels not in the other range
  gr <- suppressWarnings(c(gr1, gr2))
  gr <- subset(gr, width > min_width)
  gr <- gr[order(gr$height, decreasing=TRUE), ]
  message("Found ", length(gr), " peaks")
  gr
}
