#' Tibble of top genes from linear fit model with limma
#'
#' @param fit results from \code{lmFit}
#' @param coef contrast number
#' @param genes_columns column names in \code{fit$genes}, default is Gene.Symbol
#' and Gene.Title
#' @param res result for write_top
#' @param cutoff cutoff for write_top
#' @param \dots other options passed to \code{top_tibble}
#'
#' @return A list of tibbles
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' }
#' @export

top_tibble <- function(fit, coef=1, genes_columns){
   if(missing(genes_columns)){
      # Gene.symbol for plot_volcano  (some tables have both GENE and GENE_SYMBOL)
      if("GENE_SYMBOL" %in% colnames(fit$genes)){
        colnames(fit$genes)[colnames(fit$genes)=="GENE_SYMBOL"] <- "Gene.Symbol"
      }
      if("Gene Symbol" %in% colnames(fit$genes)){
        colnames(fit$genes)[colnames(fit$genes)=="Gene Symbol"] <- "Gene.Symbol"
      }
      ## Add SPOT_ID  - sometimes ID is not probe id!
       genes_columns <- which(colnames(fit$genes ) %in% c("Gene.Symbol",
         "Gene Title", "Gene.Title",  "GENE_NAME", "DESCRIPTION"))
   }
   if(length(genes_columns) == 0){
       # message("No columns with gene symbol or title found")
       res <- limma::topTable(fit, coef, adjust="fdr", number=Inf)
       res <- tibble::as_tibble(
                tibble::rownames_to_column(as.data.frame(res), "gene_name"))
       res <- arrange(res, gene_name)
   }else{
       res  <- limma::topTable(fit, coef, adjust="fdr", number=Inf,
                 genelist = fit$genes[ , genes_columns, drop=FALSE] )
       res <- tibble::as_tibble(
                tibble::rownames_to_column(as.data.frame(res), "probe"))
       # replace /// with ,  (or drop for title)
       if("Gene.Symbol" %in% names(res)){
           res$Gene.Symbol <-gsub(" /// ", ", ", res$Gene.Symbol)
       }
       if("Gene.Title" %in% names(res)){
           res$Gene.Title <- gsub(" //.*", "", res$Gene.Title)
       }
       res <- arrange(res, probe)
   }
   res
}


#' @describeIn top_tibble All top tibbles
#' @export
top_tibble_all <- function(fit, ...){
   c1 <- colnames( fit$contrasts)
   c1 <- gsub(" - ", "_vs_", c1)
   n <- length(c1)
   res <- vector("list", n)
   names(res) <- c1
   for(i in 1:n) res[[i]] <- top_tibble(fit, i, ...)
   res
}

#' @describeIn top_tibble Write tables to file
#' @export
write_top <- function(res, cutoff = 0.05){
     names(res) <- gsub("/", "", names(res))
     for (i in 1:length(res)) {
         vs <- gsub("\\.* ", "_", names(res[i]))
         vs <- gsub("_+_", "_", vs, fixed = TRUE)
         vs <- paste0(vs, ".txt")
        res[[i]] <- filter(res[[i]], adj.P.Val < cutoff)
        message("Saving ", nrow(res[[i]]), " rows to ", vs )
         readr::write_tsv(res[[i]], vs)
     }
}
