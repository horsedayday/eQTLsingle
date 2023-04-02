#' To visualize eQTL effect
#'
#' This function can generate violin plots of gene expression in both REF group and ALT group
#' @param eqtl.result a dataframe (only one row), describes eQTL result of one SNV-gene pair, can be generated from the function eQTLsingle
#' @param expressionMatrix dataframe of gene expression
#' @param snv.gene.pair.metadata snv-gene pair metadata, generated from the function eQTLsingle_build_metadata
#' @param figure.output optional parameter, path of saving the figure
#' @importFrom stringr str_split
#' @importFrom rlang .data
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr %>%
#' @return violin plots of gene expression in both ALT group and REF group
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
#' # Discover eQTLs
#' eQTL.result <- eQTLsingle(toy_expressionMatrix, snv.gene.pair.metadata)
#' # to draw violin plots of SNV (with snvid, "SNV_3") and gene (with geneid, "gene_16")
#' eQTL.result.target <- eQTL.result[(eQTL.result$SNVid == "SNV_3") &
#'                                   (eQTL.result$Geneid == "gene_16"),]
#' eQTLsingle_visualization(eQTL.result.target,
#'                          toy_expressionMatrix,
#'                          snv.gene.pair.metadata)
#' @export
eQTLsingle_visualization <- function(eqtl.result, expressionMatrix, snv.gene.pair.metadata, figure.output=NULL){
  SNVid <- eqtl.result$SNVid
  geneid <- eqtl.result$Geneid
  adjust.pvalue <- eqtl.result$adjusted_pvalue
  ref_cells <- snv.gene.pair.metadata[snv.gene.pair.metadata$SNVid==SNVid, 'Ref_cells']
  ref_cells <- stringr::str_split(ref_cells, ",")[[1]]
  alt_cells <- snv.gene.pair.metadata[snv.gene.pair.metadata$SNVid==SNVid, 'Alt_cells']
  alt_cells <- stringr::str_split(alt_cells, ",")[[1]]
  counts_Ref <- unlist(expressionMatrix[geneid, ref_cells])
  counts_Alt <- unlist(expressionMatrix[geneid, alt_cells])

  # build data.frame
  df <- data.frame(Expression = c(counts_Ref,counts_Alt), SNV = c(rep("REF",length(counts_Ref)), rep("ALT",length(counts_Alt))))
  violinplot <- df %>% ggplot2::ggplot(aes(color=.data$SNV, y=.data$Expression,x=.data$SNV)) +
    ggplot2::geom_violin(trim=TRUE, aes(fill = .data$SNV),color = NA) +
    ggplot2::geom_boxplot(width=0.03,outlier.shape=NA,color='#000000',lwd=0.8) +
    ggplot2::scale_fill_manual(values=c('#fd8d3c', '#807dba')) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme(axis.text=element_text(size=14)) +
    ggplot2::theme(axis.title.x = element_blank()) +
    ggplot2::theme(axis.title=element_text(size=14)) +
    ggplot2::theme(plot.margin = unit(c(1,1,1,1), "cm")) +
    ggplot2::theme(axis.title.y=element_text(vjust=4)) +
    ggplot2::theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
    ggplot2::ggtitle(paste('SNV: ', SNVid, '    ', 'Gene: ', geneid, '    ', sep = "")) +
    ggplot2::theme(plot.title = element_text(hjust=0.5,size=18))

  if(!is.null(figure.output)){
    ggplot2::ggsave(figure.output)
  }
    return(violinplot)
}
