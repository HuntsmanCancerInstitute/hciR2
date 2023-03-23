#' Plot QC Concentration vs QC RIN
#'
#' @param samples A GNomEx sample table with \code{QC RIN} and \code{QC Conc. (ng/uL)}
#' @param \dots additional options passed to hc_chart
#'
#' @return A highchart
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' plot_QC(samples, x="QC Conc. (ng/uL)", y="QC RIN")
#' }
#' @export

plot_QC <- function(samples, ...){
   highcharter::hchart(samples, "scatter", highcharter::hcaes(x = `QC Conc. (ng/uL)`, y= `QC RIN`,
              group = trt, value = ID ) ) %>%
     highcharter::hc_tooltip( pointFormat = "{point.value}", headerFormat = "") %>%
      highcharter::hc_xAxis( gridLineWidth = 1, tickLength = 0, startOnTick = "true", endOnTick = "true") %>%
      highcharter::hc_chart(zoomType = "xy", ...) %>%
      highcharter::hc_exporting(enabled=TRUE, filename = "QC")
}
