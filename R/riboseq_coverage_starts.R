#' Estimate coverage of riboseq read start around start and stop codon
#'
#' Finds coverage of the first base in the alignment between -50 and 50 bp
#' around the start codon and between -50 and 50 bp around stop codon.
#'
#' @param cds_by_tx A GRangesList from \code{cdsBy(txdb, by="tx")} with first
#' and last exon extended 50 bases using \code{resizeFeature} from
#' \code{systemPipeR}
#' @param aln GAlignments object from \code{readGAlignments}
#' @param average return entire coverage vector as RleList instead of average
#' coverage around start and stop codon
#' @param name sample name column value
#'
#' @return A tibble with start, mean start coverage (avg_start), number of
#' genes with at least one read (n_start), frame, end, avg_end, sample
#' name and number of transcripts
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  library(GenomicAlignments)
#'  aln1 <- readGAlignments(BamFile("uniq_x1.bam"))
#'  library(GenomicFeatures)
#'  txdb <- makeTxDbFromEnsembl(organism="Mus musculus", release=102)
#'  seqlevels(txdb) <- seqlevels(aln1)
#'  cds_by_tx <- cdsBy(txdb, by="tx", use.names=TRUE)
#'  cds_by_tx50 <- systemPipeR:::.resizeFeature(cds_by_tx, up=50, down=50)
#'  a1 <- subset(aln1, width(aln1) >= 27 & width(aln1) <= 30 )
#'  x1 <- riboseq_coverage_starts( cds_by_tx50, a1)
#'  x1
#' }
#' @export

riboseq_coverage_starts <- function(cds_by_tx, aln, average= TRUE,
   name = "Ribo-seq"){
   if( class(aln)[1] != "GAlignments"){
      stop("aln should be a Genomic alignment")
   }
   ## get first base of alignment only
   gr <- GenomicRanges::resize(methods::as(aln, "GRanges"), 1, fix="start")
   # calculating coverage by transcript
   cvg <- GenomicFeatures::coverageByTranscript(gr, cds_by_tx)
   ## drop transcripts with zero coverage or 1 read
   n <- sum(cvg)
   message("Dropping ", sum(n==0), " transcripts with 0 reads")
   cvg <- cvg[n > 1]
   ## drop short transcripts to avoid error (or fill matrix with zeros?)
   n <- lengths(cvg)
   if(sum(n <= 150) !=0) message("Dropping ", sum(n <= 150), " short transcripts <= 150 bp")
   cvg <- cvg[n > 150]
   if(average){
      ##  keep transcripts with at least 1 read at start ?
      cvg1 <- IRanges::heads(cvg, 151)
      n1 <- sum(cvg1)
   #   message("Dropping ", sum(n1 <2),
   #                         " transcripts with <= 1 read start near TSS")
   #   cvg <- cvg[n1>1]
   #   cvg1 <- cvg1[n1>1]
      cvg2 <- IRanges::tails(cvg, 151)
      n0 <- sum(cvg)
      n1 <- sum(cvg1)
      n2 <- sum(cvg2)

   #   sum3 <- bind_rows(complete=summary(n0), `first 100 bp` = summary(n1),
   #            `last 100 bp`=summary(n2), .id="transcript")

      message("Read starts from ", length(n0), " transcripts")
     if(length(n0) >0){

   #   print( data.frame(sum3, check.names=FALSE))
      nX1 <- -50:100
      nX2 <- -100:50

      ## AVERAGE coverage by position
      ## any conversion to vector using lapply is very slow...
      x1 <-  methods::as( methods::as(cvg1, "RleViews"), "matrix")
      x2 <-  methods::as( methods::as(cvg2, "RleViews"), "matrix")
      # how to group transcripts by gene?
      # https://plastid.readthedocs.io/en/latest/generated/plastid.bin.metagene.html
      x <- tibble::tibble(
            start = nX1,
            avg_start=colMeans(x1),
            n_start=colSums(x1 != 0),
            frame = as.factor(nX1 %% 3),
            end = nX2,
            avg_end=colMeans(x2),
            n_end=colSums(x2 != 0),
            sample=name,
            transcripts=length(n0))

        }else{
            x <- tibble::tibble()
        }
     }else{
       x <- cvg
    }
    x
}
