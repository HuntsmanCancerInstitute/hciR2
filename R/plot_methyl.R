#' Plot smoothed methylation
#'
#' Plot KEGG pathways
#'
#' @param bsseq bsseq results
#' @param dmr annotated dmrFinder output with 1 row (with start, end, and chr columns)
#' @param smooth smoothed methylation values
#' @param subset subset samples to plot
#' @param group column name for grouping, default trt
#' @param type track type
#' @param extend bases to extend up or downstream
#' @param \dots additional options passed to plotTracks
#'
#' @author Chris Stubben
#' @examples
#' \dontrun{
#'   plot_methyl
#'  }


plot_methyl <- function(bsseq, dmr, smooth=smoothPs, subset, group="trt", type =c("a", "confint"), extend=10000, ...){
   if(nrow(dmr) != 1) stop("dmr should be 1 row from annotated dmrFinder output")
   d1 <- dmr
   gtrack <- GenomeAxisTrack()
   ## title - check if columns present
   mtitle <- paste0(d1$rank, ". ", d1$name, ", ", d1$description, " (n=", d1$n,
        ", diff=", round( d1$meanDiff,3), ", stat=", round(d1$areaStat,1), ")")
   ## dataTrack
   start1 <- d1$start - extend
   end1 <- d1$end + extend
   n <- which(seqnames(gr1) == d1$chr & start(gr1) > start1 & end(gr1) < end1 )
   x1 <- gr1[n,]
   if(!length(x1) > 1) stop("no matches in gr1")

   smoothPs2 <- as.matrix(smooth[n,])
   x2 <- as.data.frame(pData(bsseq) )
   if(!group %in% colnames(x2)) stop("No column name matching ", group, ". Add the group option")

   if(!missing(subset)){
      if(!is.numeric(subset)) stop("should be a vector with column index")
       smoothPs2 <- smoothPs2[, subset]
       x2<- x2[subset,]
   }
   # factor using order in table (if subsetting to keep colors in order??)
   x2$trt <- x2[[group]]
   x2$trt <- factor(x2$trt, levels = unique(x2$trt))
   dTr <- DataTrack(data = t(smoothPs2), start = start(x1),
       end = end(x1), chromosome = seqnames(x1), name = "Methylation", genome="mm10",
       groups = x2$trt, col=unique(x2$col), legend=FALSE, type = type)

  ht <- HighlightTrack(trackList = dTr, start = d1$start, end= d1$end,
    chromosome=d1$chr, fill="NA", col="lightgrey")
  plotTracks(list(gtrack, ht, trTx), from = start1, to = end1, chromosome = d1$chr,
          transcriptAnnotation = "symbol", collapseTranscripts = "longest",
		   main = mtitle, cex.main =1, ...)
}

## plot by gene
plot_methyl2 <- function(bsseq, dmr, gene, subset, biomart=mouse96, bases=40000, shift=0, ...){
   d1 <- filter(biomart, gene_name==gene)
   if(nrow(d1)!=1) stop("No match to ", gene)
   gtrack <- GenomeAxisTrack()
   ## dataTrack
   pos <- round(mean(c(d1$start, d1$end)))
   # IF UCSC
   chr <- paste0("chr", d1$chromosome)
   # chr <-  d1$chromosome
   start1 <- pos - bases/2
   end1 <- pos + bases/2
   if(shift!=0){
      start1 <- start1 + shift
      end1 <- end1 + shift
   }
   n <- which(seqnames(gr1)==chr & start(gr1) > start1 & end(gr1) < end1 )
   x1 <- gr1[n,]

   if(!length(x1) > 1) stop("no matches in gr1")

   smoothPs2 <- as.matrix(smoothPs[n,])
   x2 <- as.data.frame(pData(bsseq) )
   if(!missing(subset)){
      if(!is.numeric(subset)) stop("should be a vector with column index")
       smoothPs2 <- smoothPs2[, subset]
       x2<- x2[subset,]
   }
   x2$trt <- factor(x2$trt, levels = unique(x2$trt))
   dTr <- DataTrack(data = t(smoothPs2), start = start(x1),
       end = end(x1), chromosome = seqnames(x1), name = "Methylation", genome="hg38",
       groups = x2$trt, col=unique(x2$col), legend=FALSE)

   d2 <- filter(dmr, name==gene)

  mtitle <- paste0(gene, " (", bases/1000, "kb window)")

  ht <- HighlightTrack(trackList = dTr, start = d2$start, end= d2$end,
    chromosome=chr, fill="NA", col="lightgrey")
  plotTracks(list(gtrack, ht, trTx), from = start1, to = end1, chromosome = chr,
          transcriptAnnotation = "symbol", collapseTranscripts = "longest",
          main = mtitle, cex.main =1, type = c("a", "confint"), ...)
}
