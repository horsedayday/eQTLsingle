#' Generate metadata for single-cell eQTL analysis.
#'
#' This function generates a dataframe describing metadata for further eQTL analysis. The metadata includes valid SNV-gene pairs and corresponding valid cells for eQTL analysis.
#' @param snvMatrix A dataframe, describes genotype of SNVs (loci by cells), rownames of this dataframe are SNVid, colnames of this dataframe are Cell ids.
#' @param expressionMatrix A dataframe, describes gene expressions (gene expressions by cells), rownames of this dataframe are Geneid, colnames of this dataframe are CellId
#' @param snv.number.of.cells threshold of minimal number of cells where mutation (or non-mutation) occurs
#' @param expression.min threshold of expression levels of valid genes, used along with another parameter, expression.number.of.cells
#' @param expression.number.of.cells threshold of minimal number of cells in which valid genes express.
#' @return A dataframe, each row describes information of valid SNV-gene pairs and corresponding valid cells for eQTL analysis
#' \itemize{
#' \item SNVid: Id of SNV, represented as CHR__POS
#' \item Num_cells_ref: number of cells in REF group
#' \item Num_cells_alt: number of cells in ALT group
#' \item Ref_cells: cell list of REF group
#' \item Alt_cells: cell list of ALT group
#' \item GeneList: gene list, these genes will be tested with the SNV in this row
#' \item Num_gene: number of genes in gene list
#' \item CellList: whole cell list of REF and ALT group
#' }
#' @examples
#' # Generate metadata for further eQTL analysis.
#' # Load snvMatrix and expressionMatrix
#' data(toy_snvMatrix)
#' data(toy_expressionMatrix)
#' # Only test the SNVs in which cells with either genotypes (REF and ALT) are at least 30 cells
#' # Only test the genes which express in > 30 cells in group REF and ALT with expression level >1
#' snv.gene.pair.metadata <- eQTLsingle_build_metadata(toy_snvMatrix,
#'                                                     toy_expressionMatrix,
#'                                                     snv.number.of.cells=30,
#'                                                     expression.min=1,
#'                                                     expression.number.of.cells=30)
#' @export
eQTLsingle_build_metadata <- function(snvMatrix, expressionMatrix, snv.number.of.cells=30, expression.min=1, expression.number.of.cells=30){
  # test if two matrix are consistence
  snvMatrix.cells <- colnames(snvMatrix)
  expressionMatrix.cells <- colnames(expressionMatrix)
  if(!setequal(snvMatrix.cells, expressionMatrix.cells)){
    stop('SNV Matrix and Gene Expression Matrix have different cells!')
  }

  snv.list <- rownames(snvMatrix)
  gene.list <- rownames(expressionMatrix)
  cell.list <- colnames(snvMatrix)

  N <- length(snv.list)
  snv.gene.pair.metadata <- data.frame(SNVid = rep("", N),
                                       Num_cells_ref = rep(0, N),
                                       Num_cells_alt = rep(0, N),
                                       Ref_cells = rep("", N),
                                       Alt_cells = rep("", N),
                                       GeneList = rep("", N),
                                       Num_gene = rep(0, N),
                                       CellList = rep("", N),
                                       stringsAsFactors = FALSE)
  useful_snv <- 0 # count snv number for final dataframe
  # snv level
  for (snvid in snv.list){
    cell.ref <- colnames(snvMatrix[snvid,snvMatrix[snvid,] == 0])
    cell.alt <- colnames(snvMatrix[snvid,snvMatrix[snvid,] == 1])
    if((length(cell.ref) > snv.number.of.cells) & (length(cell.alt) > snv.number.of.cells)){ # snv pass test
      # test gene further
      genelist <- c() # for saving genes for this snv
      for (gene in gene.list){
        cell.ref.gene <- expressionMatrix[gene, cell.ref]
        cell.alt.gene <- expressionMatrix[gene, cell.alt]
        cell.ref.gene.valid.num = sum(cell.ref.gene > expression.min) # number of valid cells on this gene
        cell.alt.gene.valid.num = sum(cell.alt.gene > expression.min) # number of valid cells on this gene
        if ((cell.ref.gene.valid.num > expression.number.of.cells) & (cell.alt.gene.valid.num > expression.number.of.cells)){ # this gene pass the test
          genelist <- c(genelist, gene)
        }
      }
      if (length(genelist) > 0){ # have valid gene
        useful_snv = useful_snv + 1 # count valid snv number
        snv.gene.pair.metadata[useful_snv, 'SNVid'] <- snvid
        snv.gene.pair.metadata[useful_snv, 'Num_cells_ref'] <- length(cell.ref)
        snv.gene.pair.metadata[useful_snv, 'Num_cells_alt'] <- length(cell.alt)
        snv.gene.pair.metadata[useful_snv, 'Ref_cells'] <- paste(cell.ref, collapse = ',')
        snv.gene.pair.metadata[useful_snv, 'Alt_cells'] <- paste(cell.alt, collapse = ',')
        snv.gene.pair.metadata[useful_snv, 'GeneList'] <- paste(genelist, collapse = ',')
        snv.gene.pair.metadata[useful_snv, 'Num_gene'] <- length(genelist)
        snv.gene.pair.metadata[useful_snv, 'CellList'] <- paste(c(cell.ref,cell.alt), collapse = ',')
      }
    }
  }
  # remove nonsense rows
  snv.gene.pair.metadata <- snv.gene.pair.metadata[1:useful_snv, ]
  return(snv.gene.pair.metadata)
}
