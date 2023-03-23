#' Add annotations to slice_with_height output
#'
#' Finds nearest feature using distanceToNearest and adds gene and biotype.
#' The distance to region will be zero if peak is within the feature.
#'
#' @param gr1 slice_to_height output
#' @param biomart biomart annotations
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'   rbga <- readGAlignments("14415X1.bam")
#'   slices <- slice_with_height(rgba, lower = 1000)
#'    slice_annotate(slices, human102)
#' }
#' @export

slice_annotate <- function(gr1, biomart){
   gr2 <- GRanges(seqnames = biomart$chromosome,
      IRanges(start = biomart$start, end = biomart$end),
      strand = biomart$strand , gene=biomart$gene_name)
   x <- IRanges::distanceToNearest(gr1, gr2)
   y <- data.frame(gr1)[, c(1,2,5,4,7)]
   z <- cbind(y, biomart[subjectHits(x), 2:3])
   colnames(z)[1]<- "chr"
   ## distance to range?
   z$to_range <- data.frame(x)$distance
   ## or distance to start
   y1 <- start(gr1[S4Vectors::queryHits(x), ])
   y2 <- start(gr2[S4Vectors::subjectHits(x), ])
   z$to_start  <- y1 - y2
   as_tibble(z)
}
