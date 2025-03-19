#' Cluster -log10 pvalues from top gene sets in each contrast
#'
#' @param e1 results from enrichr_all
#' @param db List index or name, default 1
#' @param top_n Number of top gene sets to plot
#' @param trim Trim long names, default 50
#' @param min_p Minimum p-value to plot (to adjust scale)
#' @param complexHeatmap use pheatmap from ComplexHeatmap
#' @param \dots additional options passed to pheatmap
#'
#' @return pheatmp
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' library(enrichR)
#' setEnrichrSite("Enrichr")
#' celldbs <- c("KEGG_2021_Human", "WikiPathways_2019_Human")
#' e1 <- enrichr_all(res, celldbs)
#' plot_top_enrichr(e1, 1, top=3)
#' }
#' @export

plot_top_enrichr <- function(e1, db=1, top_n=5, trim=50, min_p=100, complexHeatmap=TRUE, ...){
  if(class(e1)[1] == "list"){
     if(grepl("^[0-9]+$", db)) db <- names(e1)[db]
     if(is.na(db)) stop("No databases found, use a name or number < ", length(e1)+1)
     e1 <- e1[[db]]
     if(is.null(e1)) stop("No match to ", db)
     # change cluster to contrast?
     if("cluster" %in% names(e1)) e1 <- dplyr::rename(e1, contrast="cluster")
      message("Plotting top ", top_n, " sets from ", db)
  }
  z <- dplyr::filter(e1, Adjusted.P.value == 0 )
  if(nrow(z) > 0){
    minP <-  min(e1$Adjusted.P.value[e1$Adjusted.P.value !=0])
    message("Warning: ", paste(unique(z$Term), collapse=", "),
     " with Adjusted p-value = 0, setting to ", scales::scientific(minP, digits = 1) )
    e1$Adjusted.P.value[e1$Adjusted.P.value == 0] <- minP
  }
  y <- dplyr::group_by(e1, contrast) %>%
		     dplyr::slice_min(Adjusted.P.value, n=top_n, with_ties=FALSE)
  # get all clusters with top n terms
	y <- dplyr::filter(e1, Term %in% y$Term)
  y$Term <- gsub(" \\(GO:.*", "", y$Term)
  y$Term <- gsub("[ _]WP[0-9]+$", "", y$Term)
  z <- dplyr::select(y, contrast, Term, Adjusted.P.value) %>%
	      dplyr::mutate(Adjusted.P.value = -log10(Adjusted.P.value)) %>%
	    	tidyr::spread(contrast, Adjusted.P.value)
  clrs <- grDevices::colorRampPalette(
			 c("white", RColorBrewer::brewer.pal(n = 9, name = "Reds")[-1]))(255)
  ## too many NAs to cluster
  z <- as_matrix(z)
  z[is.na(z)] <- 0
  # fix -log10(0) = Inf ...
  z[is.infinite(z)] <- max(z[!is.infinite(z)])
  if(any(z > min_p)){
	  message("Setting minimum p-value to 1e-", min_p)
	z[z > min_p] <- min_p
  }
  # wrapping in pheatmap creates a large gap
  rownames(z) <- ifelse(nchar(rownames(z)) > trim,
				 paste0(substr(rownames(z), 1, trim-2), "..."), rownames(z))
  if(complexHeatmap){
    ComplexHeatmap::pheatmap(z,  color = clrs,
      heatmap_legend_param = list(title = "-log10\npvalue", tick_length = unit(0, "mm"),
      legend_width = unit(50, "mm")) , ...)
  }else{
    pheatmap::pheatmap(z, color = clrs, ...)
  }
}
