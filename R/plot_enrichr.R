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
#' @param cluster_row cluster rows
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

plot_enrichr <- function(e1, db=1, vs=1, plot="bar", width=0.7, top=10, wrap=40,
   xlab_bar="pathway", top_n, cluster_row=FALSE, exclude_row, exclude_col, set_first=FALSE, ...){
   if(class(e1)[1] != "list") stop("enrichr_all results should be a list")
   if(grepl("^[0-9]+$", db)) db <- names(e1)[db]
   if(is.na(db)) stop("No databases found, use a name or number < ", length(e1)+1)
   e1 <- e1[[db]]
   if(is.null(e1)) stop("No match to ", db)
   # cluster?
   if("cluster" %in% names(e1)) e1 <- dplyr::rename(e1, contrast="cluster")
   if(!is.character(vs)) vs <- unique(e1$contrast)[vs]
   if(is.na(vs)) stop("No contrast found, use a name or number < ", length(unique(e1$contrast))+1 )
   if(! vs %in% unique(e1$contrast)) stop("No contrast matching ", vs)
   ## format Terms
   e1$Term <- gsub(" \\(GO:.*", "", e1$Term)
   e1$Term <- gsub("[ _]WP[0-9]+$", "", e1$Term)
   e1$Term <- wrap_string(e1$Term, wrap)
   #-------------------------------------------------------
   # if topN then plot all contrasts in image grid
   if(!missing(top_n)){
	    message("Plotting top ", top_n, " sets in each contrast from ", db)
	    z <- dplyr::group_by(e1, contrast) %>%
		       dplyr::slice_min(Adjusted.P.value, n=top_n, with_ties=FALSE)
	    if(!is.factor(z$Term)) z$Term <- factor(z$Term, levels = rev(sort(unique(z$Term))))
      ## fix 0 p-values
         minP <-  min(z$Adjusted.P.value)
         if(minP==0){ message("Note: setting p-values = 0 to 1e-100 for log10 scale")
               minP <- 1e-100
               z$Adjusted.P.value[z$Adjusted.P.value == 0] <- 1e-100
         }
         ggplot2::ggplot(z, aes(x=contrast, y=Term, fill=Adjusted.P.value)) +
         ggplot2::geom_tile(width=0.9, height=0.9) +
         ggplot2::xlab("") +
         ggplot2::ylab("") +
         ggplot2::theme_classic() +
         ggplot2::theme(
            axis.text.x = ggplot2::element_text(size=11, angle = 45, hjust=1),
            axis.text.y = ggplot2::element_text(size=11),
            legend.key =  ggplot2::element_rect(size=1.5, colour="white")) +
         ggplot2::scale_fill_gradient(low="#A50F15", high="#FEE0D2", trans = "log10", na.value="white",
             limits=c(minP, max(z$Adjusted.P.value))) +
         ggplot2::guides(fill=ggplot2::guide_legend(title="Adjusted\n p-value", reverse = FALSE))
   }else{
      message("Using contrast ", vs, " from ", db)
      y <- dplyr::filter(e1, contrast == vs) %>%
           dplyr::arrange(Adjusted.P.value)
      if(!missing(exclude_row)) y <- dplyr::filter(y, !Term %in% exclude_row)
      if(nrow(y) == 0) stop("No significant results found")
      y <- utils::head(y, top)  %>% dplyr::arrange(desc(row_number()))
      #-------------------------------------------------------
      ## Cluster genes
      if(tolower(plot) %in% c("cluster", "genes")){
         x <- dplyr::select(y, Term, Genes)  %>%
	             tidyr::separate_rows(Genes, sep=";") %>%
		           dplyr::mutate(Term = forcats::fct_rev(forcats::fct_inorder(Term)))  %>%
               dplyr::mutate(present=1)
      x <- tidyr::spread(x, Genes, present, fill=0)
      callback <- function(hc, ...) {  dendsort::dendsort(hc) }
	    # cluster_row = FALSE by default??
      pheatmap::pheatmap(as_matrix(x), col=c("white","red"), legend=FALSE, border=NA,
            cluster_rows=cluster_row, clustering_callback = callback, ...)
   }else{
      y <- dplyr::mutate(y, N=as.integer(gsub("/.*", "", Overlap))) %>%
	           dplyr::arrange(desc(Adjusted.P.value), N) %>%
             dplyr::mutate(Term = forcats::fct_inorder(Term))
     #-------------------------------------------------------
 	   ## GEOM_TILE contrasts
     if(tolower(plot) %in% c("tile", "contrasts")){
		    z <- dplyr::filter(e1, Term %in% y$Term)
        if(!missing(exclude_col)) z <- dplyr::filter(z, !contrast %in% exclude_col)
		    z$Term <- factor(z$Term, levels = y$Term)
		    ## set first contrast
		    if(!is.factor(z$contrast)) z$contrast <- factor(z$contrast)
		    if(set_first) z$contrast <- forcats::fct_relevel(z$contrast, vs)
        ggplot2::ggplot(z, aes(x=contrast, y=Term, fill=Adjusted.P.value)) +
          ggplot2::geom_tile(width=0.9, height=0.9) +
          ggplot2::xlab("") +
          ggplot2::ylab("") +
          ggplot2::theme_classic() +
          ggplot2::theme(
             axis.text.x = ggplot2::element_text(size=11, angle = 45, hjust=1),
             axis.text.y = ggplot2::element_text(size=11),
             legend.key = ggplot2::element_rect(size=1.5, colour="white")) +
          ggplot2::scale_fill_gradient(low="#A50F15", high="#FEE0D2", trans = "log10", na.value="white",
             limits=c(min(z$Adjusted.P.value), max(z$Adjusted.P.value))) +
          ggplot2::guides(fill=ggplot2::guide_legend(title="Adjusted\n p-value", reverse = FALSE))
      }else{
        #-------------------------------------------------------
        # Bar
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
    }
  }
}
