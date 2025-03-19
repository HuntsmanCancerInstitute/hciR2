#' Loop through FindAllMarkers and run EnrichR
#'
#' @param clust results from FindAllMarkers with cluster and gene columns
#' @param celldbs a vector of library names from listEnrichrDbs
#' @param padj Adjusted P value cutoff for Enrichr results
#'
#' @return a list of tibbles
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' library(enrichR)
#' setEnrichrSite("Enrichr")
#' celldbs <- c("KEGG_2021_Human", "WikiPathways_2019_Human")
#' res <- FindAllMarkers(object = srt, only.pos = TRUE, logfc.threshold = 0.5)
#' e1 <- enrichr_clusters(res, celldbs)
#' }
#' @export

enrichr_clusters <- function(clust, celldbs, padj = 0.1){
   x1 <- unique(clust$cluster)
   n1 <- length(x1)
   c1 <- vector("list", n1)
   names(c1) <- x1
   for (i in 1:n1){
      message("Check cluster ", x1[i])
      sig <- unique(sort(clust$gene[clust$cluster == x1[i]]))
      e1 <- enrichr(sig, celldbs)
      c1[[i]] <- e1
   }
   n2 <- length(celldbs)
   e1 <- vector("list", n2)
   names(e1) <- celldbs
   for(i in 1:n2){
	  # no results returns data.frame with all logical columns
	  z <- purrr::map(c1, celldbs[i])
    z2 <- bind_rows(z[sapply(z, nrow) > 0], .id = "cluster") %>% as_tibble()
     if(nrow(z2) == 0){
         message("WARNING: no matches to ", celldbs[i], ". Please check spelling")
       } else{
         z2 <-  filter(z2, Adjusted.P.value < padj)
       }
       z2$cluster <- factor(as.integer(z2$cluster))
      e1[[i]] <- z2
   }
   e1
}
