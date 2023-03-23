#' Plot multiple enrichment plots
#'
#' @param x a table from fgsea_all (top n rows are displayed)
#' @param fc fold change vector
#' @param db gene set database
#' @param n number of rows in x to display, default 4
#' @param ncol number of columns in plot_gird, default 2
#' @param contrast add contrast name to title
#' @param scaleplots scale option for plot_grid, default 0.95
#' @param sameyaxis use the same y-axis range for all plots, default TRUE
#'
#' @return a plot_grid
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' multi_enrich_plot(h1[[1]], fc[[1]], msig_hallmark, n=6)
#' }
#' @export

multi_enrich_plot  <- function(x, fc, db, n=4, contrast, ncol=2,  scaleplots = 0.95, sameyaxis=TRUE){
  ## format subtitle with NES and Padj
  x$NES <- format(round(x$NES, digits=2), nsmall = 2)
  x$p <- ifelse(x$padj < 0.01, formatC(x$padj, format = "e", digits = 0), formatC(x$padj, digits = 2, flag="#") )
  p <- vector("list", n)
  for(i in 1:n){
    # use NULL if missing
    if(is.na(x$pathway[i])){
    # R FAQ - you can set elements to NULL using x[i] <- list(NULL)
	  p[i] <- list(NULL)
	}else{
    t1 <- x$pathway[i]
    if(!missing(contrast)) t1 <- paste(contrast, t1)
    p[[i]] <- plotEnrichment(db[[x$pathway[i] ]], fc) +
    labs(title=t1, subtitle = paste0("NES=", x$NES[i], ", Padj=", x$p[i]) ) +
    xlab("Rank") + ylab("Enrichment score")
	}
  }
  if(sameyaxis){
     # apply this new range to all plot objects
	 nn <- !sapply(p, is.null)
     n2 <- range(sapply(p[nn], function(x) layer_scales(x)$y$range$range))
     p[nn] <- lapply(p[nn], function(x) x + ylim(n2[1], n2[2] ))
  }
  cowplot::plot_grid(plotlist=p, ncol=ncol, scale = scaleplots) + theme(plot.margin=margin(10,10,10,10))
}


# plot the same pathway in all contrasts.  May need to run fgsea without a cutoff

# k1 <- fgsea_all(res, msig_pathways$KEGG, FDR=1)
# z <- bind_rows(k1, .id = "contrast") %>% filter(pathway=="Apoptosis")

multi_enrich_contrast  <- function(x, fc, db, n=4, ncol=2,  scaleplots = 0.95, sameyaxis=TRUE){
  ## format subtitle with NES and Padj
  x$NES <- format(round(x$NES, digits=2), nsmall = 2)
  x$p <- ifelse(x$padj < 0.01, formatC(x$padj, format = "e", digits = 0), formatC(x$padj, digits = 2, flag="#") )
  p <- vector("list", n)
  for(i in 1:n){
    # use NULL if missing
    if(is.na(x$pathway[i])){
    # R FAQ - you can set elements to NULL using x[i] <- list(NULL)
	  p[i] <- list(NULL)
	}else{
    t1 <- x$contrast[i]
    p[[i]] <- plotEnrichment(db[[x$pathway[i] ]], fc[[t1]]) +
    labs(title=t1, subtitle = paste0("NES=", x$NES[i], ", Padj=", x$p[i]) ) +
    xlab("Rank") + ylab("Enrichment score")
	}
  }
  if(sameyaxis){
     # apply this new range to all plot objects
	 nn <- !sapply(p, is.null)
     n2 <- range(sapply(p[nn], function(x) layer_scales(x)$y$range$range))
     p[nn] <- lapply(p[nn], function(x) x + ylim(n2[1], n2[2] ))
  }
  cowplot::plot_grid(plotlist=p, ncol=ncol, scale = scaleplots) + theme(plot.margin=margin(10,10,10,10))
}




# PVALUE
multi_enrich_plot_pvalue  <- function(x, fc, db, n=4, ncol=2, scaleplots = 0.95, sameyaxis=TRUE){
  ## format subtitle with NES and Padj
  x$NES <- format(round(x$NES, digits=2), nsmall = 2)
  x$p <- ifelse(x$pval < 0.01, formatC(x$pval, format = "e", digits = 1), formatC(x$pval, digits = 2, flag="#") )
  p <- vector("list", n)
  for(i in 1:n){
    # use NULL if missing
    if(is.na(x$pathway[i])){
    # R FAQ - you can set elements to NULL using x[i] <- list(NULL)
	p[i] <- list(NULL)
	}else{
    p[[i]] <- plotEnrichment(db[[x$pathway[i] ]], fc) +
    labs(title=x$pathway[i], subtitle = paste0("NES=", x$NES[i], ", P-value=", x$p[i]) ) +
    xlab("Rank") + ylab("Enrichment score")
	}
  }
  if(sameyaxis){
     # apply this new range to all plot objects
	 nn <- !sapply(p, is.null)
     n2 <- range(sapply(p[nn], function(x) layer_scales(x)$y$range$range))
     p[nn] <- lapply(p[nn], function(x) x + ylim(n2[1], n2[2] ))
  }
  cowplot::plot_grid(plotlist=p, ncol=ncol, scale = scaleplots) + theme(plot.margin=margin(5,5,5,5))
}

# option to plot more than 4?  set by=n from plot_multi_enrich
gsea4pdf <- function(h1, fc, db, filename, pos=8, neg=8){
 x1 <- filter(h1, enriched=="positive") %>% head(pos)
 x2 <- filter(h1, enriched=="negative") %>% head(neg)
 pdf(file=filename, height=8, width=10.5)
  for(i in seq(1,pos, by=4)) print(multi_enrich_plot(x1[i:(i+4),], fc, db))
  for(i in seq(1,neg, by=4)) print(multi_enrich_plot(x2[i:(i+4),], fc, db))
 dev.off()
}
