#' Run enricher from clusterProfiler on findAllMarkers
#'
#' @param res A list of DESeq results
#' @param gsets A list of gene sets (or two column dataframe with pathway and gene)
#' @param maxGSSize Maximum set size, default 2000
#' @param padj Adjusted p-value cutoff for significant genes, default 0.05
#' @param quiet suppress messages except total significant
#' @param \dots Additional options like minGSSize and pvalueCutoff passed to \code{enricher}
#'
#' @return A list of tibbles
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' library(hciRdata)
#' enricher_all(res, msig_pathways$KEGG)
#' }
#' @export

enricher_markers <- function(res, gsets, maxGSSize = 2000, padj = 0.05, quiet=TRUE, ...){
   ## convert to 2 column data.frame
   if(class(gsets)[1] == "list"){
      term2gene <- dplyr::bind_rows(
		 lapply(gsets, tibble::enframe, value="SYMBOL"), .id="pathway") %>%
         dplyr::select(1,3)
  }else{
	  term2gene <- gsets
  }
   if(!is.data.frame(res)) stop("res should be a data.frame with gene and cluster columns from findAllMarkers")
   res <- filter(res, p_val_adj < padj)
   # hack to fix gene_name
   n1 <- which(names(res)=="gene_name")
   if(length(n1)==1)  names(res)[n1] <- "gene"
   res$cluster <- as.character(res$cluster)

   clusters <- unique(res$cluster)
   n <- length(clusters)
   x <- vector("list", n)
   names(x) <-  clusters

   for(i in 1:n){
       r1 <- filter(res, cluster == clusters[i])
	   ## enricher will return error if no genes in set!
       if(any(r1$gene %in% term2gene$SYMBOL)){

   	     em1 <-  clusterProfiler::enricher(r1$gene, TERM2GENE=term2gene, maxGSSize=maxGSSize, ...)
	     if(nrow(em1) > 0){
	        # ID and Description are the same
	        z <- tibble::as_tibble(em1)[, -2]
	        names(z)[1] <- "pathway"
            x[[i]]  <- z
	    }
		if(!quiet)	message("Cluster ", clusters[i], ", ", nrow(z), " enriched sets")
	 }
   }
   x1 <- bind_rows(x, .id="cluster")
   x1$cluster <- factor(x1$cluster, levels=unique(x1$cluster))
   x1
}
