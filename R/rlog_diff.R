#' Calculate fold changes using rlog differences
#'
#' @param rld a DESeqTransform object with regularized log values
#' @param exp column numbers for treatment
#' @param ref column numbers for reference level, will be repeated if needed
#' @param biomart annotations from read_biomart
#'
#' @return a list of result tables with log2FoldChange
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'   rlog_diff(rld, "G1", mouse96)
#' }
#' @export


rlog_diff <- function(rld, exp, ref, biomart){
   if(class(rld) != "DESeqTransform") stop("rld should be a DESeqTransform object")
   if(length(ref) < length(exp) ) ref <- rep(ref, length=length(exp))
   x <- SummarizedExperiment::assay(rld)
   res <- vector("list", length(exp))
   names(res) <- paste(colnames(x)[exp], colnames(x)[ref], sep = " vs. ")
   ids <- rownames(x)
   for(i in seq_along(exp)){
      fc <-  x[, exp[i]] - x[, ref[i]]
      ## include biotype?
      a1 <- biomart[ ,c(1,2,4,8)]
      if("human_homolog" %in% names(biomart)) a1 <- biomart[ ,c(1,2,4,8,10)]
      y <- inner_join( a1,
        tibble(id = ids, x[, exp[i]], x[, ref[i]], log2FoldChange = fc), by="id") %>%
        arrange(desc(log2FoldChange))
      names(y)[ncol(y) - 2:1 ] <- colnames(x)[c(exp[i], ref[i])]
       res[[i]] <- y
  }
  res
}
