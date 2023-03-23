#' Estimate riboseq coverage for different read sizes.
#'
#' Loops through an alignment file and runs code{riboseq_coverage_start}
#'
#' @param lengths a vector of alignment lengths
#' @param cds A GRangesList from \code{cdsBy(txdb, by="tx")} with first and last exon
#' extended 50 bases using \code{resizeFeature} from \code{systemPipeR}
#' @param aln GAlignments object from \code{readGAlignments}
#'
#' @return A tibble with read length and columns from code{riboseq_coverage_start}
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  psite1 <- psite_offsets(26:31, cds_by_tx50, aln1)
#'  dplyr::select(psite1, 1,8) %>% unique()
#'  ggplot(psite1, aes(x=start, y=avg_start, color=length)) + geom_line() +
#'  xlab("Start of read relative to start codon") +
#'  ylab("Reads per transcript") + theme(legend.position="none") +
#'  facet_wrap( ~ length, ncol=1, scales="free", strip.position="right")
#' }
#' @export

riboseq_psite_offsets <- function(lengths, cds, aln){
   n <- length(lengths)
   x <- vector("list", n)
   names(x) <- lengths
   for(i in seq_along(lengths)){
      message("Coverage for read length ", lengths[i])
	  # use njunc = 0?
      aX <- subset(aln,  width == lengths[i] )
      #print(length(aX))
      x[[i]] <- suppressMessages( riboseq_coverage_starts( cds, aX))
  }
   bind_rows(x, .id="length")
}
