#' List samples in Expression sets
#'
#' @param eset expression set from \code{GEOquery}
#' @param unique show columns with same repeated values
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @export

eset_samples <- function(eset, unique=FALSE){
  if(class(eset) != "ExpressionSet") stop( "eset should be an ExpressionSet")
  p1 <- Biobase::pData(eset)
  if(unique){
     n1 <- apply(p1, 2, function(x) length(unique(x)))
     message("Listing ", sum(n1==1), " columns with the exact same values")
     x <- t(p1[1, n1==1])
     p1 <- tibble( rowname = rownames(x), value =as.vector(x[,1]) )
  }else{
  ## RENAME characteristics_ch1* columns
  n <- grep("characteristics_ch1", names(p1))
   if(length(n) > 0){
      ## OR column ending in :ch1 , see GSE46268
      n1 <- grep(":ch1$", names(p1))
      if(length(n1) == length(n)){
         message("Dropping :ch1 from ", paste(names(p1)[n1], collapse=", "))
         names(p1)[n1] <- gsub(":ch1", "", names(p1)[n1])
         p1 <- p1[, -n]
      }else{
        for(i in n){
           p1[,i] <- as.character(p1[,i])
           new_col_name <- unique(gsub( " ", "_", gsub(":.*", "", p1[,i])))
           new_col_name <- new_col_name[new_col_name!=""]
           if(length(new_col_name) != 1){
              message("Note: cannot parse name in mixed column ", i)
              next
           }
           names(p1)[i] <- new_col_name
           p1[,i] <- gsub(".*: ", "", p1[,i])
       }
     }
  }
  ## drop columns with same value
  n1 <- apply(p1, 2, function(x) length(unique(x)))
  message("Dropping ", sum(n1==1), " columns with the exact same values
  (add unique=TRUE to view these columns)")
  p1 <- p1[, n1 > 1]
  # drop ftp supplementary_file?
  n <- grepl("supplementary_file", names(p1))
  if(any(n)){
     message("Dropping supplementary_file columns with long ftp strings")
     p1 <- p1[ , !n ]
  }
  ## add row.name - same as geo_accession
  n2 <- sapply(p1, function(x) identical(x, rownames(p1)))
  if(any(n2)){
     message(" Note: row names are same as ", names(which(n2)))
  }else{
     message("Adding row names to column 1")
     p1 <- cbind( rowname = rownames(p1), p1)
  }
  p1 <- lapply(p1, as.character)
  p1 <- tibble::as_tibble(p1)
  p1 <- suppressMessages( readr::type_convert(p1))
  }
  p1
}
