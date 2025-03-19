#' Plot KEGG pathway maps using pathview
#'
#' Plot KEGG pathways
#'
#' @param x a results table with gene_name and log2FoldChange
#' @param pathway_name match pathway name
#'
#' @author Chris Stubben
#' @examples
#' \dontrun{
#'   data(kegg_hsa)
#'   data(egSymb)   # in gage
#'   plot_kegg(res, "Cytokine")
#'  }
#' @export

plot_kegg <- function(x, pathway_name){
 y <- unique(kegg_hsa[, 1:3]) %>% filter(grepl(pathway_name, pathway))
 if(nrow(y)==0) stop("No match to ", pathway)
 if(nrow(y)> 1) message("Matches ", nrow(y), "pathways, using first")
 message("Plotting ", y$pathway[1])
 keggid <- y$entry[1]
 # x is annotated results with gene_name and log2FoldChange
 # convert to named vector with sorted fold changes
  y <- suppressMessages( write_gsea_rnk(x, write=FALSE))
  rnk <- y$log2FoldChange
  names(rnk) <- y$gene_name
  names(rnk) <- gage::sym2eg(names(rnk))
   pathview::pathview(rnk, pathway.id = keggid)
}
