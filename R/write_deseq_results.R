#' Write DESeq results
#'
#' Write DESeq result files to an Excel file.
#'
#' @param result_all a list from \code{results_all}
#' @param sets table with set intersections, optional
#' @param file file name
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  write_deseq_results(res)
#' }
#' @export

write_deseq_results <- function(result_all, sets, file = "DESeq.xlsx"){

   ##  if results are a tibble (since simplify=TRUE by default)
   if(!class(result_all)[1] == "list"){
         n <- attr(result_all, "contrast")
         result_all <- list(result_all)
         names(result_all) <- n
   }
   ## 1. summary
   sum1 <-  dplyr::bind_rows(lapply(result_all, summary_deseq), .id= "contrast")
   # write.xlsx does not like tibbles, so use as.data.frame
   sum1 <- as.data.frame(sum1)
   # write.xlsx requires a named list for writing mulitple worksheets
   sum1 <- list("summary" = sum1)
   res1 <- lapply(result_all, as.data.frame )

   ## write.xlsx replaces space with .
   names(res1)  <- gsub( "\\.? ", "_", names(res1))
   ## forward slash will cause Excel errors
   names(res1)  <- gsub( "/", "", names(res1))

   DESeq_tables <-  c( sum1, res1)
   if(!missing(sets)) DESeq_tables <-  c( sum1, res1, list(Sets=sets))
   message("Saving ", length(DESeq_tables), " worksheets to ", file)
  # DESeq_tables
   openxlsx::write.xlsx(DESeq_tables, file = file, rowNames= FALSE )

}
