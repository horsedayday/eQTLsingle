#' Generate SNV matrix
#'
#' This function generates SNV matrix (loci by cells) from the output of the function snv_readcounts_tidy. In this SNV matrix, (1) -1 indicates missing data (no sufficient reads covered); (2) 0 indicates genotype of this locus is REF (non-mutated); (3) 1 indicates genotype of this locus is ALT (mutated);
#' @param readcountProfile A dataframe describes number of reads of each base (A,C,G,T) over the locus. This dataframe can be generated from the function snv_readcounts_tidy
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
#' @return A dataframe, describes genotypes of each locus (loci by cells), (1) -1 indicates missing data (no sufficient reads covered); (2) 0 indicates genotype of this locus is REF (non-mutated); (3) 1 indicates genotype of this locus is ALT (mutated)
#' @export

snv_matrix <- function(readcountProfile){
  cell.whole <- unique(unlist(readcountProfile[,"CELL"]))
  readcountProfile$SNVid = paste(readcountProfile$CHROM, readcountProfile$POS, sep = "__")
  snvid.whole <- unique(unlist(readcountProfile[,"SNVid"]))
  # build a data.frame snvid * cell
  snvMatrix <- data.frame(matrix(-1, nrow = length(snvid.whole), ncol = length(cell.whole)))
  rownames(snvMatrix) <- snvid.whole
  colnames(snvMatrix) <- cell.whole
  for (j in 1:dim(readcountProfile)[1]){ # each row
    snvid <- readcountProfile[j, 'SNVid']
    cell <- readcountProfile[j, 'CELL']
    genotype <- readcountProfile[j, 'Genotype']
    snvMatrix[snvid, cell] <- ifelse(genotype == 'REF', 0, 1) # set genotype value, ref-0, alt-1
  }
  return(snvMatrix)
}
