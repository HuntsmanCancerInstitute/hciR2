#' Barplot of top pathways, regulators, disease or toxicity functions.
#'
#' Bars are colored by z-score and height is p-value.  If stacked=TRUE, plot the
#' percent of up and down-regulated genes or activated and inhibited regulators by Category
#'
#' @param ipa a table from read_IPA
#' @param top Number of top pathways, regulators or functions to display, default 15
#' @param pvalue P-value cutoff, default 0.05
#' @param wrap wrap pathway strings, default 40 characters
#' @param barwidth geom_bar width, default 0.7
#' @param name_category add Category to x-axis name, default FALSE
#' @param stacked plot stacked barchart, default FALSE
#' @param results DESeq results loaded into IPA, required for stacked pathway barplot
#'
#' @return a ggplot
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' ipa1 <- read_ipa("IPA_export_all.txt")
#' barplot_ipa(ipa1[[2]])
#' barplot_ipa(ipa1[[2]], stacked=TRUE, results=res[[2]])
#' barplot_ipa(ipa1[[3]], name_cat=TRUE, wrap=30)
#' }
#' @export

barplot_ipa <- function(ipa, top = 15, pvalue = 0.05, wrap = 30, barwidth = 0.7,
   name_category =FALSE, stacked = FALSE, results){
   p1 <- check_ipa(ipa)
   # default IPA pathway barplot
   if (!stacked) {
      y <- head(p1, top)  %>%
	       arrange(desc(row_number()))
      if(name_category) y$Pathway <- paste0(y$Pathway, " (", y$Categories, ")" )
      y <- mutate(y, Pathway = wrap_string(Pathway, wrap)) %>%
        mutate(Pathway = fct_inorder(Pathway))
      n1 <- max(abs(y$zScore), na.rm=TRUE)
      ggplot(y, aes(x = Pathway, y = -log10(pValue), fill = zScore)) +
   	    geom_bar(width = barwidth, stat = "identity", color = "gray70") +
   	    xlab("") +
		    ylab("-log10(P-value)") +
   	    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
   	    theme_bw() +
		    theme(axis.text.y=element_text(size=11)) +
		    coord_flip() +
   	    scale_fill_gradient2(name = "Z-score", low = "blue", mid = "white", high = "orange", limits=c(-n1, n1))
   } else {
     ## STACKED barplot - pathways
      if(!"Categories" %in% names(p1)){
        if(missing(results)) stop("Stacked barchart with pathways requires DESeq results with log2FoldChange")
        if("human_homolog" %in% names(results)) results$gene_name <- results$human_homolog
        # remove duplicates
        results <- arrange(results, padj) %>% filter(!duplicated(gene_name))
        z <- head(p1, top)  %>%
           arrange(desc(row_number())) %>%
             mutate(Pathway = wrap_string(Pathway, wrap)) %>%
             mutate(Pathway = fct_inorder(Pathway)) %>%
             dplyr::select(Pathway, Molecules)
        z2 <- separate_rows(z, Molecules, sep=",") %>%
        	inner_join(results, by=c(Molecules = "gene_name" ))%>%
        	group_by(Pathway) %>% summarize(n = n(), up = sum(log2FoldChange > 0)/n, down = 1-up) %>%
          mutate(Pathway = factor(Pathway, levels = z$Pathway)) %>%
          dplyr::select(Pathway, n, up, down)
        # extra space for total counts
        ymax <- 1.09
        maxn <- max(z2$n)
        if(maxn > 99) ymax <- 1.12
        z3 <- gather(z2, "regulated", "total", -(1:2))
        p1 <- ggplot(z3, aes(x=Pathway, y=total, fill=regulated, group=regulated)) +
          geom_bar(width=barwidth, stat="identity", color="gray70",
                  position = position_stack(reverse = TRUE)) +
          xlab("") +
          ylab("Significant genes (%)") +
          scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0))) +
          theme_classic() +
          theme(axis.text.y=element_text(size=11)) +
          coord_flip()  +
          scale_fill_manual(values=c("green", "red"))
        p1 +  geom_text(data=subset(z3, regulated=="up"), aes(x=Pathway, y=1.01, label=n),
            size=4, color="black", hjust="left")
      } else{
       ## STACKED barplot - regulators and diseases
      # Separate comma-separated list of Categores into new rows
      # and Dont split DNA Replication, Recombination, and Repair
      xlab <- "Regulators"
      if (names(p1)[1] == "Categories") {
		     xlab <- "Diseases"
         p1 <- mutate(p1, Categories = gsub(", Recombination,", " Recombination", Categories)) %>%
               separate_rows(Categories , sep=", ")
      }
	    #option to filter by abs z-score > 1
      y <-  filter(p1, pValue < pvalue, abs(zScore)>0) %>%
	           group_by(Categories)  %>%
			       mutate(Categories = wrap_string(Categories, wrap)) %>%
             summarize(  n = n(),
                         activated = sum(zScore > 0, na.rm=TRUE)/n ,
                         inhibited =  1 - activated  ) %>%
	           arrange(desc(n))
      # extra space for total counts
      ymax <- 1.09
      maxn <- max(y$n)
      if(maxn > 99) ymax <- 1.12
      if(maxn > 199) ymax <- 1.14
      z <- head(y, top) %>% gather("state", "total", -(1:2))  %>%
             mutate(Categories = factor(Categories, levels=y$Categories[top:1]))
      p1 <- ggplot(z, aes(x=Categories, y=total, fill=state, group=state)) +
        geom_bar(width=0.7, stat="identity", color="gray70",
                position = position_stack(reverse = TRUE)) +
	     xlab("") +
	     ylab(xlab) +
	     scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0))) +
	     theme_classic() +
	     theme(legend.title = element_blank(), axis.text.y=element_text(size=11)) +
       coord_flip()  +
       scale_fill_manual(values=c("orange", "blue"))
      p1 +  geom_text(data=subset(z, state=="inhibited"), aes(x=Categories, y=1.01, label=n),
           size=4, color="black", hjust="left")

      }
   }
}
