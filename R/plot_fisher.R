#' Plot Fisher all output
#'
#' Plot heatmap with p-values or percent by contrast
#'
#' @param x list from \code{link{fgsea_all}}
#' @param trim trim long names, default more than 70 characters
#' @param sets display contrasts sharing n or more sets for n > 1.  If n = 1,
#' then only plot unique sets.  If missing, then plots all sets, default.
#' @param p.value pvalue p-value cutoff
#' @param percent plot percent in set, default is -log10(pvalue)
#' @param palette color palette, default Reds
#' @param cluster_row Cluster dendrogram rows, default is an alphabetical list
#' @param \dots other options passed to \code{pheatmap}
#' @author Chris Stubben
#' @examples
#' \dontrun{
#'   library(hciRdata)
#'   x <- fisher_all(fc, msig_pathways$KEGG)
#'   plot_fisher(x)
#' }
#' @export

plot_fisher <- function(x, trim=70, sets, p.value=0.05, percent=FALSE, palette="Reds", cluster_row=FALSE, ...){
   if(is.data.frame(x)) stop("A list from fisher_all is required")
   y <- dplyr::bind_rows(x, .id = "contrast") %>%
       filter(pvalue < p.value)
   ## order columns by order in list (or alphabetical)
   y$contrast <- factor(y$contrast, levels= names(x))
   y$term <- ifelse(nchar(y$term) > trim,
                   paste0(substr(y$term, 1, trim-2), "..."), y$term)
   ## pvalue or percent
   if(!percent){
   z <- dplyr::select(y, contrast, term, pvalue) %>%
        dplyr::mutate(pvalue = -log10(pvalue)) %>%
         tidyr::spread(contrast, pvalue)
   }else{
   z <- dplyr::select(y, contrast, term, percent) %>%
         tidyr::spread(contrast, percent)
   }
   if(!missing(sets)){
     n <- apply(z[, -1], 1, function(x) sum(!is.na(x)))
     if(sets ==1){
        z <- filter(z, n == 1)
     }else{
        z <- filter(z, n >= sets)
     }
   }
   clrs <- palette255(palette)
   ## too many NAs to cluster
   z <- as_matrix(z)
   z[is.na(z)] <- 0
   message(nrow(z) , " total sets")
   n1 <- max(abs(z), na.rm=TRUE)
   brks <- seq(0, n1, length = 255)
   pheatmap::pheatmap(z, color = clrs, breaks = brks, cluster_cols=FALSE,
       cluster_rows=cluster_row, ...)
}
