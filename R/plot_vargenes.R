#' Bar plot of variance fractions for top genes
#'
#' @param x result from fitExtractVarPartModel in variancePartition package
#' @param column column name
#' @param n number of top variable genes, default 20
#' @param barwidth barplot pointWidth
#' @param first put column first in key
#' @param \dots additional options passed to \code{hc_size}
#'
#' @return A highchart
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  x <- read_tsv("biotypes.txt")
#'  plot_biotypes(x)
#' }
#' @export

plot_vargenes <- function(x, column, n = 20, barwidth=NULL, first=FALSE, ...){
   if(!column %in% names(x)) stop("Column is missing")
   if(is_tibble(x)){
      # gene name in column 2?
      x$gene_name[x$gene_name==""] <- x$id[x$gene_name==""] 
      x <- as_matrix( x[,-1])
   }
   x <- as.data.frame(x)
   colNames <- colnames(x)
   n1 <- order(x[[column]], decreasing=TRUE)
   y <- head(x[n1,], n)
   y <- round(y*100, 2)
   z <- tibble::rownames_to_column(y, "gene") %>%
     gather( "col", "var", -gene)
   if(first){
     # column first and residuals last
     n2 <- sort(unique(z$col))
     z$col <- factor(z$col, levels = c(column, n2[!n2 %in% c( column, "Residuals")], "Residuals"))
   }else{
      z$col <- factor(z$col, levels = colNames )
   }
   highcharter::hchart(z, "bar", highcharter::hcaes(x=gene, y= var, group = col)) %>%
   highcharter::hc_plotOptions(series = list(stacking = stack),
      bar=list(pointWidth=barwidth)) %>%
   highcharter::hc_yAxis(title=list(text="Variance explained (%)"), max=100, reversedStacks=FALSE) %>%
   highcharter::hc_xAxis(title=list(text="")) %>%
   highcharter::hc_plotOptions( series=list(states=list(inactive=list(opacity=1)))) %>%
   highcharter::hc_size(...) %>%
   highcharter::hc_exporting(enabled=TRUE, filename = "vargenes")
}
