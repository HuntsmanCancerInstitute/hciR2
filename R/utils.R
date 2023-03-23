# fix Long reactome names
cut_string <- function(x, length = 40){
   x <- trimws(x)
   n <- nchar(x) > length
   x <- substr(x, 1, length)
   x[n] <- gsub(" [^ ]*$", "", x[n])
   x <- gsub(",$", "", x)
   x
}

# wrap long names
wrap_string <- function(x, length = 50){
   paste0(substr(x, 1, length-1), sub(" ", "\n", substring(x, length)))
}

# p-value dot size using https://rdrr.io/github/surh/HMVAR/src/R/plotgg_manhattan.r
minus_log10 <- scales::trans_new("minus_log10",
                  transform = function(x) -log10(x),
                  inverse   = function(x) 10^(-x),
                  breaks    = scales::log_breaks(base = 10),
                  domain    = c(0,Inf))

# check IPA table
check_ipa <- function(ipa){
  if (class(ipa)[1] == "list" ){
	 message("NOTE: ipa is a list, using pathways")
	 p1 <- ipa[["Canonical Pathways"]]
  } else{
	 p1 <- ipa
   }
  if(nrow(p1)==0) stop("No rows found")
  c1 <- names(p1)[1]
  if(c1 == "Upstream Regulator"){
	  names(p1)[c(1,3,5,9,10)] <- c("Pathway", "Categories", "zScore", "pValue", "Molecules")
	  # N only needed for barplots if fill = p-value
	  p1$N <- stringr::str_count(p1$Molecules, ",") + 1
  }else if(c1 == "Master Regulator"){
	   names(p1)[c(1,2,7,8,11)] <- c("Pathway", "Categories", "zScore", "pValue", "Molecules")
	   p1$N <- stringr::str_count(p1$Molecules, ",") + 1
  }else if(c1 == "Categories"){
	  names(p1)[c(3,4,6)] <- c("Pathway", "pValue", "zScore")
	  p1$N <- p1$`# Molecules`
  }else if(c1 != "Pathway"){
	  stop("Cannot match first column name to pathway, regulator or function tables")
  }
  p1
}
