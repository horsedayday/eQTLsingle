#' Tidy results of eQTLsingle-SNVPostAnalysis.sh
#'
#' This function is to tidy output from the script eQTLsingle-SNVPostAnalysis.sh into a data frame.
#'
#' @param readcounts.result The output of eQTLsingle-SNVPostAnalysis.sh
#' @param minimal.read.num The threshold of minimal number of reads of each base (A, C, G, T) over the locus
#' @return
#' A data frame, describes number of reads of each base over different loci in single cells
#' \itemize{
#' \item CHROM: chromosome
#' \item POS: position
#' \item REF: reference base
#' \item ALT: alternative base
#' \item DP: depth (number of reads over the locus)
#' \item ADepth: number of reads with base A over the locus
#' \item CDepth: number of reads with base C over the locus
#' \item GDepth: number of reads with base G over the locus
#' \item TDepth: number of reads with base T over the locus
#' \item CELL: cell
#' \item Genotype: a binary variable, describing the genotype of the locus. 'REF' indicates the genotype of the locus is non-mutated (all of reads over it are same as reference). 'ALT' indicates the genotype of this locus is mutated (some of reads over it are same as alternative allele)
#' }
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom tidyr separate
#' @importFrom magrittr %>%
#' @importFrom utils read.table
#' @export

snv_readcounts_tidy <- function(readcounts.result, minimal.read.num = 5){
  readcountProfile <- utils::read.table(readcounts.result,
                                 stringsAsFactors = FALSE, col.names=c('CHROM', 'POS', 'REF', 'DP',
                                                                       'NON1', 'Ainfo', 'Cinfo', 'Ginfo', 'Tinfo',
                                                                       'NON2', 'Other', 'CELL'), fill=TRUE)
  readcountProfile$NON1 <- NULL
  if.cell.is.blank <- readcountProfile$CELL == ""
  readcountProfile$CELL[if.cell.is.blank] <- readcountProfile$Other[if.cell.is.blank]
  readcountProfile$Other[if.cell.is.blank] <- ""
  readcountProfile <- readcountProfile %>% tidyr::separate(.data$Ainfo, c(NA,'ADepth',rep(NA, 12)), sep=":")
  readcountProfile <- readcountProfile %>% tidyr::separate(.data$Cinfo, c(NA,'CDepth',rep(NA, 12)), sep=":")
  readcountProfile <- readcountProfile %>% tidyr::separate(.data$Ginfo, c(NA,'GDepth',rep(NA, 12)), sep=":")
  readcountProfile <- readcountProfile %>% tidyr::separate(.data$Tinfo, c(NA,'TDepth',rep(NA, 12)), sep=":")
  readcountProfile <- readcountProfile %>% tidyr::separate(.data$NON2, c(NA,'NDepth',rep(NA, 12)), sep=":")
  options(warn=-1)
  readcountProfile <- readcountProfile %>% tidyr::separate(.data$Other, c('OtherBase','OtherDepth',rep(NA, 12)), sep=":")
  options(warn=1)
  readcountProfile$ADepth <- as.numeric(readcountProfile$ADepth)
  readcountProfile$CDepth <- as.numeric(readcountProfile$CDepth)
  readcountProfile$GDepth <- as.numeric(readcountProfile$GDepth)
  readcountProfile$TDepth <- as.numeric(readcountProfile$TDepth)
  readcountProfile$NDepth <- as.numeric(readcountProfile$NDepth)
  readcountProfile$OtherDepth <- as.numeric(readcountProfile$OtherDepth)
  readcountProfile$OtherDepth[is.na(readcountProfile$OtherDepth)] <- 0
  readcountProfile$DP <- readcountProfile$DP - readcountProfile$NDepth - readcountProfile$OtherDepth
  readcountProfile <- readcountProfile %>% dplyr::filter(.data$NDepth == 0 & .data$OtherDepth==0)
  readcountProfile$OtherBase <- NULL
  readcountProfile$OtherDepth <- NULL
  readcountProfile$NDepth <- NULL
  readcountProfile <- readcountProfile %>% dplyr::filter((.data$ADepth >= minimal.read.num| .data$ADepth ==0) &
                                                    (.data$CDepth >= minimal.read.num | .data$CDepth ==0) &
                                                    (.data$GDepth >= minimal.read.num | .data$GDepth ==0) &
                                                    (.data$TDepth >= minimal.read.num | .data$TDepth ==0))
  readcountProfile$REF <- toupper(readcountProfile$REF)
  readcountProfile <- readcountProfile %>%
    dplyr::mutate(Genotype = case_when(.data$REF == 'A' ~ case_when(.data$ADepth == .data$DP ~ 'REF',
                                                       TRUE ~ 'ALT'),
                                       .data$REF == 'G' ~ case_when(.data$GDepth == .data$DP ~ 'REF',
                                                       TRUE ~ 'ALT'),
                                       .data$REF == 'T' ~ case_when(.data$TDepth == .data$DP ~ 'REF',
                                                       TRUE ~ 'ALT'),
                                       .data$REF == 'C' ~ case_when(.data$CDepth == .data$DP ~ 'REF',
                                                       TRUE ~ 'ALT')))
  options(warn=-1)
  snv.ref_alt.reference <- readcountProfile %>%
    dplyr::group_by(.data$CHROM, .data$POS) %>%
    dplyr::summarise_each(funs(sum),-.data$CELL,-.data$REF,-.data$Genotype) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(readcountProfile %>%
                        dplyr::select(.data$CHROM, .data$POS, .data$REF) %>%
                        unique(), by=c('CHROM', 'POS'))
  options(warn=1)

  snv.ref_alt.reference$ALT <- NA # add ALT information
  snv.ref_alt.reference$ALT <- as.character(snv.ref_alt.reference$ALT)

  for (i in 1:dim(snv.ref_alt.reference)[1]){
    ref <- snv.ref_alt.reference[i,'REF']
    refCol <- paste0(ref, 'Depth')
    df.tmp <- unlist(snv.ref_alt.reference[i, c('ADepth', 'CDepth', 'GDepth', 'TDepth')])
    order.tmp <- order(-df.tmp)[1:2]
    ref_and_alt <- names(df.tmp)[order.tmp]
    altCol <- ref_and_alt[ref_and_alt!=refCol]
    alt <- unlist(strsplit(altCol,""))[1]
    snv.ref_alt.reference[i,'ALT'] <- alt
  }

  readcountProfile <- inner_join(readcountProfile, snv.ref_alt.reference %>% dplyr::select(.data$CHROM, .data$POS, .data$REF, .data$ALT), by=c('CHROM', 'POS','REF'))
  rm(snv.ref_alt.reference)
  readcountProfile.colnames <- colnames(readcountProfile)
  readcountProfile.colnames <- c(readcountProfile.colnames[1:3],'ALT',readcountProfile.colnames[4:10])
  readcountProfile <- readcountProfile[,readcountProfile.colnames]
  return(readcountProfile)
}
