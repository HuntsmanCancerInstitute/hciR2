#' Rename biomart table to SAF format for featureCounts
#'
#' @param x biomart table
#' @param gene_name use gene_name as GeneID, default Ensembl ID
#'
#' @return tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' }
#' @export

biomart_to_SAF <- function(x, gene_name=FALSE){
  if(gene_name){
     x$GeneID <- ifelse(x$gene_name=="", x$id, x$gene_name)
     n <- sum(duplicated(x$GeneID))
     message(n, " duplicated gene names, fixing with make.unique")
     x$GeneID <- make.unique(x$GeneID)
  }else{
     x$GeneID <- x$id
  }
   rename(x, chromosome="Chr", start="Start", end="End", strand="Strand") %>%
   mutate(Strand = ifelse(Strand==1, "+", "-")) %>%
   dplyr::select(GeneID, Chr, Start,End,Strand)
}
