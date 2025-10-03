#' #' Run Enrichr on all DESeq2 result tables
#'
#' @param res A list of DESeq results
#' @param celldbs A vector of gene set names from listEnrichrDbs
#' @param deseq.padj cutoff for significant genes to test, default 0.05
#' @param fisher.padj cutoff for testing overlap, default 0.05
#' @param logFC log2FoldChange cutoff for significant genes, default none
#' @param protein_coding Test protein coding genes only, default TRUE
#'
#' @note if Enrichr is not responding, check SSL certificate
#'  httr::set_config(httr::config(ssl_verifypeer = 0L))
#'
#' @return tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' library(enrichR)
#' setEnrichrSite("Enrichr")
#' celldbs <- c("KEGG_2021_Human", "WikiPathways_2019_Human")
#' e1 <- enrichr_all(res, celldbs)
#' }
#' @export

enrichr_all <- function(res, celldbs, deseq.padj = 0.05, fisher.padj = 0.05, logFC,
  protein_coding = TRUE){
  ##  if results are a tibble (since simplify=TRUE by default)
  if(!inherits(res, "list")){
		 n <- attr(res, "contrast")
		 if(is.null(n)) stop("Missing contrast name")
		 res <- list(res)
		 names(res) <- n
  }
  n <- length(res)
  enrich <- vector("list", n)
  vs <- names(res)
  names(enrich) <-  vs
  for(i in 1:n){
	   message( i, ". ", vs[[i]])
	   s1 <- dplyr::filter(res[[i]], padj <= deseq.padj, !is.na(gene_name))
     n1a <- nrow(s1)
	   if(protein_coding){
        message("Dropping non-coding genes")
        s1 <- dplyr::filter(s1, biotype == "protein_coding")
     }
     if(!missing(logFC)) s1 <- dplyr::filter(s1, abs(log2FoldChange) >= logFC)
     if(nrow(s1) == 0){
		    message("  No significant genes")
		    enrich[[i]]  <- tibble::tibble()
     }else{
         sig <- unique(s1$gene_name)
         n2a <- n1a - length(sig)
		     message("  ", length(sig), " significant genes (dropped ", n2a, " duplicate or missing gene names)")
         # print enrichR output for first result table only
		     if(i == 1){
             x <- enrichR::enrichr(sig, celldbs)
         }else{
             ## suppress cat
             invisible(capture.output( x <- enrichR::enrichr(sig, celldbs)))
         }
          y <- lapply(x, function(y) as_tibble(filter(y, Adjusted.P.value< fisher.padj )))
          # no results causes error with bind_rows
          y[sapply(y, nrow)==0] <- NULL
         enrich[[i]] <- y
     }
  }
  n2 <- length(celldbs)
  x  <- vector("list", n2)
  names(x) <- celldbs
  for(i in 1:n2){
      x[[i]] <- dplyr::bind_rows(purrr::map(enrich, celldbs[i]), .id="contrast")  
  }
  x
}
