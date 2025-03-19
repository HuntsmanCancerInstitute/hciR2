#' Loop through FindMarkers and run EnrichR
#'
#' @param res two or more results from FindMarkers, combined into a single table with contrast column
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
#' e1 <- enrichr_contrasts(res, celldbs)
#' }
#' @export

enrichr_contrasts <- function(res, celldbs, padj = 0.1){
   vs <- unique(res$contrast)
   n1 <- length(vs)
   c1 <- vector("list", n1)
   names(c1) <- vs
   for (i in 1:n1){
      message("Check contrast ", vs[i])
      sig <- unique(sort(res$gene[res$contrast == vs[i]] ))
      e1 <- enrichr(sig, celldbs)
      c1[[i]] <- e1
   }
   n2 <- length(celldbs)
   e1 <- vector("list", n2)
   names(e1) <- celldbs
   for(i in 1:n2){
	    # no results returns data.frame with all logical columns
	    z <- purrr::map(c1, celldbs[i])
      z2 <- bind_rows(z[sapply(z, nrow) > 0], .id = "cluster") %>%
            as_tibble()
      if (nrow(z2) == 0) {
            message("WARNING: no matches to ", celldbs[i], ". Please check spelling")
      }
      else {
        z2 <- filter(z2, Adjusted.P.value < padj)
      }
      e1[[i]] <- z2
   }
   e1
}

#' @describeIn enrichr_contrasts Run Enrichr on single FindMarkers result
#' @export
enrichr_gene <- function(res, celldbs, padj = 0.1){
      sig <- unique(sort(res$gene))
      e1 <- enrichr(sig, celldbs)
      e2 <- bind_rows(e1, .id="database") %>% as_tibble() %>%
         filter(Adjusted.P.value < padj)
      e2
}
