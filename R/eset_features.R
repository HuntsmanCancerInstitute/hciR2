#' List features in Expression sets
#'
#' @param eset expression set from \code{GEOquery}
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @export

eset_features <- function(eset){
   if(class(eset) != "ExpressionSet") stop( "eset should be an ExpressionSet")
   f1 <- Biobase::fData(eset)
   ##
   # n1 <- apply(f1, 2, function(x) length(unique(x)))
   f1 <- lapply(f1, as.character)
   f1 <- tibble::as_tibble(f1)
   f1 <- suppressMessages( readr::type_convert(f1))
   f1
}
