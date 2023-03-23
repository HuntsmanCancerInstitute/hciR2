#' Cluster genes in top IPA pathways
#'
#' @param ipa a table from read_IPA
#' @param top Number of top pathways, regulators or functions to display, default 7
#' @param cut Cut long pathway names, default 50 characters
#' @param \dots additional options passed to scale_x_continuous
#'
#' @return a ggplot
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' ipa1 <- read_ipa("IPA_export_all.txt")
#'
#' }
#' @export

cluster_ipa <- function(ipa, top=10, cut=50, ...){
   p1 <- check_ipa(ipa)
   y <- head(p1, top)
   ## wrapping strings does not work with pheatmap...
   x <- dplyr::select(y, Pathway, Molecules)  %>%
	    separate_rows(Molecules, sep=",") %>%
		mutate(Pathway = cut_string(Pathway, cut),
      Pathway = fct_inorder(Pathway))  %>%
        mutate(present=1)
   x <- spread(x, Molecules, present, fill=0)
   callback <- function(hc, ...) {  dendsort::dendsort(hc) }
   # cluster_row = FALSE by default
   pheatmap::pheatmap(as_matrix(x), col=c("white","red"), legend=FALSE, border=NA,
      clustering_callback = callback, ...)
}
