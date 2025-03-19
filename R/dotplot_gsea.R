#' Dotplot for fgsea_all output with normalized enrichment scores and dots sized
#' by adjusted p-values
#'
#' @param fgsea_all_res a tibble from fgsea_all
#' @param top Number of top pathways by NES to display, default 7
#' @param pad padding for x-axis, default 0.1
#' @param wrap wrap pathway strings, default 50 characters
#' @param \dots additional options passed to scale_x_continuous
#'
#' @return a ggplot
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' h1 <- fgsea_all(res, msig_hallmark)
#' dotplot_gsea(h1[[1]])
#' # pass labels to scale_x_continuous if too many digits
#' dotplot_gsea(h1[[1]], labels = function(x) format(round(x, digits=2), nsmall = 2))
#' }
#' @export

dotplot_gsea <- function(fgsea_all_res, top=7, pad=0.1, wrap=50, ...){
   if(class(fgsea_all_res)[1] != "tbl_df") stop("fgsea_all_res should be a table from fgsea_all, try adding an index [[1]]?")
   ## reorder since fgsea_all sorted by padj then NES..
   dot_df <- filter(fgsea_all_res, enriched == "positive") %>%
               head(top) %>% arrange(desc(NES)) %>%
               bind_rows(filter(fgsea_all_res, enriched == "negative") %>%
	           head(top) %>% arrange(desc(NES)))
  dot_df$pathway <- wrap_string(dot_df$pathway, wrap)
  mutate(dot_df, pathway = fct_inorder(pathway)) %>%
    ggplot(aes(x = NES, y = pathway, size = padj)) +
    geom_point() +
    xlab("Normalized Enrichment Score") +
    ylab("")  +
    theme_bw() +
    theme(axis.text.y=element_text(size=10)) +
    scale_size("Adjusted\n p-value", trans = minus_log10) +
    facet_grid(. ~ enriched, scales = "free") +
    scale_x_continuous(expand = expansion(mult = c(pad, pad)), ...)
}
