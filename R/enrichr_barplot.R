#' Plot enrichr_all output
#'
#' @param e1 results from enrichr_all
#' @param db List index or name, default 1
#' @param vs Contrast number or name, default 1
#' @param plot plot type, cluster, tile or default bar
#' @param width bar width
#' @param top Number of top sets to plot
#' @param wrap Wrap long string names, default 40
#' @param xlab_bar xlab name
#' @param top_n plot top sets in each contrast
#' @param cluster_rows cluster rows
#' @param exclude_row exclude rows
#' @param exclude_col exclude columns
#' @param set_first set as first
#' @param \dots additional options passed to pheatmap
#'
#' @return plots
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' library(enrichR)
#' setEnrichrSite("Enrichr")
#' celldbs <- c("KEGG_2021_Human", "WikiPathways_2019_Human")
#' e1 <- enrichr_all(res, celldbs)
#' plot_enrichr(e1, 1,2)
#' plot_enrichr(e1, 1,2, plot="cluster", fontsize_col=7)
#' }
#' @export

enrichr_barplot <- function(e1, vs=1, width=0.7, top=10, wrap=40, xlab_bar="pathway", ...){
   if(!inherits(e1, "tbl_df")) stop("e1 should be a single result table from enrichr_all")
   # wrap long names
  e1$Term <- gsub("[ _]WP[0-9]+$", "", e1$Term)
   e1$Term <- wrap_string(e1$Term, wrap)


   if(!is.character(vs)) vs <- unique(e1$contrast)[vs]
   if(is.na(vs)) stop("No contrast found, use a name or number < ", length(unique(e1$contrast))+1 )
   if(! vs %in% unique(e1$contrast)) stop("No contrast matching ", vs)

   message("Using contrast ", vs)

   y <- dplyr::filter(e1, contrast == vs) %>%
        dplyr::arrange(Adjusted.P.value)

   if(nrow(y) == 0) stop("No significant results found")
   y <- utils::head(y, top)  %>% dplyr::arrange(desc(row_number()))

   y <- dplyr::mutate(y, N=as.integer(gsub("/.*", "", Overlap))) %>%
          dplyr::arrange(desc(Adjusted.P.value), N) %>%
          dplyr::mutate(Term = forcats::fct_inorder(Term))

    ## 0 p-values return scale error  'from' must be a finite number
    minP <-  min(y$Adjusted.P.value)
    if(minP == 0){ message("Note: setting p-values = 0 to 1e-100 for log10 scale")
       minP <- 1e-100
       y$Adjusted.P.value[y$Adjusted.P.value == 0] <- 1e-100
    }
    ggplot2::ggplot(y, aes(x=Term, y=N, fill=Adjusted.P.value)) +
      ggplot2::geom_bar(width=width, stat="identity") +
      ggplot2::xlab("") +
      ggplot2::ylab(paste("Significant genes in", xlab_bar)) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.05))) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.y=element_text(size=11)) +
      ggplot2::scale_fill_gradient(low="#A50F15", high="#FEE0D2", trans = "log10",
        limits=c(minP,  max(y$Adjusted.P.value) )) +
      ggplot2::guides(fill=ggplot2::guide_legend(title="Adjusted\n p-value", reverse = FALSE))

}
