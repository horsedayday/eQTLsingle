% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snv_reference.R
\name{snv_reference}
\alias{snv_reference}
\title{Generate genotype reference table of SNVs}
\usage{
snv_reference(readcountProfile)
}
\arguments{
\item{readcountProfile}{A dataframe describes number of reads of each base (A,C,G,T) over the locus. This dataframe can be generated from the function snv_readcounts_tidy
\itemize{
\item CHROM: chromosome
\item POS: position
\item REF: reference base
\item ALT: alternative base
\item DP: depth (number of reads over the locus)
\item ADepth: number of reads with base A over the locus
\item CDepth: number of reads with base C over the locus
\item GDepth: number of reads with base G over the locus
\item TDepth: number of reads with base T over the locus
\item CELL: cell
\item Genotype: a binary variable, describing the genotype of the locus. 'REF' indicates the genotype of the locus is non-mutated (all of reads over it are same as reference). 'ALT' indicates the genotype of this locus is mutated (some of reads over it are same as alternative allele)
}}
}
\value{
A data frame, describes base types of reference allele and alternative allele of each SNV.
\itemize{
\item SNVid: location of SNV, represented as CHR__POS
\item REF: reference base
\item ALT: alternative base
}
}
\description{
This function generates a genotype reference table, i.e., base types of reference allele and alternative allele of each SNV. Input of this function can be generated from the function snv_readcounts_tidy
}
