#' Convert enrichr_all to enrichResults
#'
#' https://github.com/jmw86069/multienrichjam/blob/629b347911ae3306d8a6ae1007668ec86a9a292f/R/jamenrich-base.r
#'
#' @param e1 results from enrichr_all
#' @param db List index or name, default 1
#' @param vs Contrast number or name, default 1
#' @param genes vector of significant genes tested
#'
#' @return enrichResults object
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' library(enrichR)
#' setEnrichrSite("Enrichr")
#' celldbs <- c("KEGG_2021_Human", "WikiPathways_2019_Human")
#' e1 <- enrichr_all(res, celldbs)
#' sig <- sort(unique(res[[2]]$gene_name[res[[2]]$padj < 0.05]))
#' e2 <- as_enrichResult(e1, 1,2, genes = sig)
#' library(ClusterProfiler)
#' barplot(e2, show=10)
#' dotplot(e2)
#' heatplot(e2, show=10)
#' library(enrichplot)
#' e2 <- pairwise_termsim(e2)
#' emapplot(e2, show=10)
#' emapplot(e2, show=10, cex_line=3) + ggraph::scale_edge_width(range = c(1,10))
#' }
#' @export

as_enrichResult <- function(e1, db=1, vs=1, genes){
   if(class(e1)[1] != "list") stop("enrichr_all results should be a list")
   if(grepl("^[0-9]+$", db)) db <- names(e1)[db]
   if(is.na(db)) stop("No databases found, use a name or number < ", length(e1)+1)
   e1 <- e1[[db]]
   if(is.null(e1)) stop("No match to ", db)
   if(grepl("^[0-9]+$", vs)) vs <- unique(e1$contrast)[vs]
   if(is.na(vs)) stop("No contrast found, use a name or number < ", length(unique(e1$contrast))+1 )
   message("Formatting contrast ", vs, " from ", db, " as enrichResult")
   e1 <- filter(e1, contrast == vs) %>% arrange(Adjusted.P.value)

   ## only GO and Wikipathway have IDs at end of term???
   ids <- gsub(".*_", "", e1$Term)
   ## format Terms
   e1$Term <- gsub(" \\(GO:0.*", "", e1$Term)
   e1$Term <- gsub("_WP[0-9]+$", "", e1$Term)

   z <- data.frame(
 	  ID          = ids,
      Description = e1$Term,
      GeneRatio   = e1$Overlap,
      pvalue      = e1$P.value,
      p.adjust    = e1$Adjusted.P.value,
      geneID      = gsub(";", "/", e1$Genes),
	  Count       = as.integer(gsub("/.*", "", e1$Overlap)),
	  row.names   = ids)

   x1 <- new("enrichResult",
      result        = z,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      gene          = genes,
      universe      = "UNKNOWN",
      geneSets      = list(sets=1:5),
      organism      = "Human",
      keytype       = "UNKNOWN",
      ontology      = db,
      readable      = FALSE)
   x1
}
