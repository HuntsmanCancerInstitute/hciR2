#' Get a count matrix of top probes
#'
#' Join DESeq results to rlog or other counts, mainly for plotting gene heatmaps.
#'
#' @param res results from top_tibble
#' @param eset expression set with normalized counts
#' @param padj adjusted p-value cutoff, default 0.05
#' @param n number of top probes, default 40
#' @param row_names row names, Gene.Symbol default
#' @param by join results to count rownames by this column, default probe
#'
#' @return A tibble with colData attribute
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  x <- top_counts(res[[1]], eset)
#'  pheatmap(as_matrix(x))
#'  plot_genes(x, "trt")
#' }
#' @export

top_probes <- function(res, eset,  padj = 0.05, n=40, row_names = "Gene.Symbol", by="probe"){
   if(!class(eset) == "ExpressionSet") stop("eset shoud be a ExpressionSet")
   esetx <- Biobase::exprs(eset)
   colx <- Biobase::pData(eset)
   x <- dplyr::filter(res, adj.P.Val < padj)
    if(nrow(x)==0) stop("No matching rows with padj < ", padj)
    if(nrow(x)==1) stop("Only one matching row, nothing to plot")
   x <- dplyr::arrange(x, adj.P.Val) %>% utils::head(n)
   if( !row_names %in% colnames(x)) row_names <- "gene_name"
   ## match column 1 in results to count rownames
   n <- match(x[[ by ]], rownames(esetx))
   if(all(is.na(n))) stop("Column ", by, " in results and rownames in counts do not match")
   if(any(is.na(n))) stop(sum(is.na(n)) , " result rows not in count matrix" )
   mat <- esetx[n, ]
   mat <- dplyr::bind_cols( tibble::tibble(id=x[[ row_names ]]), tibble::as_tibble(mat) )
   attr(mat, "colData") <- colx
   mat
}
