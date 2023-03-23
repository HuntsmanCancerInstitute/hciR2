#' Read all Excel spreadsheets
#'
#' @param file Excel file name
#'
#' @return A list of tibbles
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'    read_excel_files()
#' }
#' @export

read_excel_sheets <- function(file){
   s1 <- readxl::excel_sheets(file)
   x <- lapply(s1, function(n) readxl::read_excel(file, sheet=n ))
   names(x) <- s1
   x
}
