#' Read STAR SJ.out  file
#'
#' @param file file name
#'
#' @return tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' }
#' @export

# column 1: chromosome
# column 2: first base of the intron (1-based)
# column 3: last base of the intron (1-based)
# column 4: strand (0: undefined, 1: +, 2: -)
# column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
# column 6: 0: unannotated, 1: annotated (only if splice junctions database is used)
# column 7: number of uniquely mapping reads crossing the junction
# column 8: number of multi-mapping reads crossing the junction
# column 9: maximum spliced alignment overhang

read_sjout <- function(file){
 x <- read_delim(file, col_names= c("chr", "start", "end", "strand", "motif",
     "annotation", "unique", "multimap", "overhang"))
 x
}
