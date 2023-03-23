#' Dotplot for IPA output with normalized enrichment scores and dots sized
#' by adjusted p-values
#'
#' @param ipa a table from read_IPA
#' @param top Number of top pathways, regulators or functions to display, default 7
#' @param pvalue P-value cutoff, default 0.05
#' @param wrap wrap pathway strings, default 40 characters
#' @param pad padding for x-axis, default 0.1
#' @param first plot inhibited values first (down diagonal)
#' @param \dots additional options passed to scale_x_continuous
#'
#' @return a ggplot
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' ipa1 <- read_ipa("IPA_export_all.txt")
#' dotplot_ipa(ipa1)
#' }
#' @export


dotplot_ipa <- function(ipa, top=7, pvalue = 0.05, wrap=30, pad=0.1, first="inhibited", ...){
   p1 <- check_ipa(ipa)
   # arrange significant pathways by z-score
   y <- filter(p1, !is.na(zScore), pValue < pvalue) %>%
         arrange(zScore) %>%
         mutate(state = factor(ifelse(zScore < 0, "inhibited", "activated"),
           levels=c("inhibited", "activated")))
   dot_df <- filter(y, state == "inhibited") %>%
         head(top) %>%
         bind_rows( filter(y, state == "activated") %>%
		 tail(top))
   dot_df$Pathway <- wrap_string(dot_df$Pathway, wrap)
   if (first == "inhibited") dot_df <- dot_df[nrow(dot_df):1, ]
   mutate(dot_df, Pathway = fct_inorder(Pathway)) %>%
     ggplot(aes(x=zScore, y=Pathway, size=pValue)) +
      geom_point() +
	  xlab("Activation Z-Score") +
	  ylab("")  +
      theme_bw() +
	  theme(axis.text.y=element_text(size=11)) +
      scale_size("P-value", trans=minus_log10)+
	  facet_grid(.~state, scales="free") +
      scale_x_continuous(expand = expansion(mult = c(pad, pad)), ...)
}
