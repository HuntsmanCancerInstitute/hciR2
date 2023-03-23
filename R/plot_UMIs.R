#' Plot UMIs found by qiagen_smallRNA_umi_extractor
#'
#' @param x file with sample id, UMIs found and no UMI tags
#' @param stack, plot percent or total reads (percent or normal)
#' @param width bar width, default 15
#' @param \dots additional options passed to \code{hc_size}
#'
#' @return A highchart
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  x <- read_tsv("UMIs.txt")
#'  plot_UMIs(x)
#' }
#' @export

plot_UMIs <- function(x, stack = "percent", width = 15, ...){
   names(x)[1] <- "sample"
   if(stack == "percent"){
       ylab <- "Percent"
   }
   else{
      stack <- "normal"
      ylab <- "Total Reads"
   }
   df <- tidyr::gather(x, "feature", "count", -sample)
   df$feature <- ifelse(df$feature=="UMI", "UMI found", "No UMI")
   # PLOT percent or total
   highcharter::hchart(df, "bar",
      highcharter::hcaes(x=sample, y= count, group = feature)) %>%
   highcharter::hc_plotOptions(series = list(stacking = stack),
      bar=list(pointWidth=width)) %>%
   highcharter::hc_yAxis(title=list(text=ylab), reversedStacks=FALSE) %>%
   highcharter::hc_xAxis(title=list(text="")) %>%
   highcharter::hc_size(...) %>%
   highcharter::hc_exporting(enabled=TRUE, filename = "UMIs")
}
