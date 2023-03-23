#' Plot heatmap for a single chromosome
#'
#' Plot pheatmap for genes on a single chromosome and labels every 10, 20, 30 Mb
#' and optionally add a vector of gene names
#'
#' @param rld rlog or other counts in DESeqTransform object
#' @param biomart filter annotations from \code{read_biomart}
#' @param gene_names a vector of gene names to label
#' @param \dots other options passed to \code{pheatmap}
#'
#' @return A pheatmap
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  pheatmap_chr(rld,  subset(mmu, chromosome==11),  gene_names = "Ewsr1" )
#' }
#' @export

pheatmap_chr <- function(rld, biomart, gene_names,  ...){
   x <- assay(rld)
   chr <- arrange(biomart, start)
   n <- match( chr$id, rownames(x))
   y <- x[n,]
   ## add zero expression for filtered counts
   y[is.na(y)] <- 0
   rownames(y) <- rep(" ", nrow(y))
   ## print 10 Mb, 20 Mb, etc at closest gene
   n <-table(floor(chr$start/1e7)*10)
   rownames(y)[cumsum(n[-length(n)])+1]  <-  paste(names(n[-1]), "Mb")
   ## print some gene names
   if(!missing(gene_names)){
      n <- which( chr$gene_name %in% gene_names)
      if(length(n) > 0) rownames(y)[n] <-  chr$gene_name[n]
   }
   colrs <- rev((grDevices::colorRampPalette(
                  RColorBrewer::brewer.pal(11, "RdBu")))(255))
   pheatmap(y, col=colrs, cluster_rows=FALSE,  border=NA, scale="row", ...)
}
