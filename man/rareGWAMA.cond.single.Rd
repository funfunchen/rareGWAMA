% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rareGWAMA.conditional.single.R
\name{rareGWAMA.cond.single}
\alias{rareGWAMA.cond.single}
\title{conditional analysis for single variant association test;}
\usage{
rareGWAMA.cond.single(
  score.stat.file,
  imp.qual.file = NULL,
  ref.file,
  candidateVar,
  knownVar,
  alternative = "two.sided",
  ...
)
}
\arguments{
\item{score.stat.file}{the file names of score statistic files;}

\item{imp.qual.file}{the file names of imputation quality;}

\item{ref.file}{the file names of the reference panel file; For VCF references, use the VCF file name; for binary PLINK references, use the prefix;}

\item{candidateVar}{the tabix range;}

\item{knownVar}{known variant;}

\item{alternative}{The alternative hypothesis. Default is two.sided;}

\item{col.impqual}{The column number for the imputation quality score;}

\item{impQual.lb}{The lower bound for the imputation quality. Variants with imputaiton quality less than impQual.lb will be labelled as missing;}

\item{impQualWeight}{Using imputation quality as weight}

\item{rmMultiAllelicSite}{Default is TRUE. Multi-allelic sites will be removed from the analyses if set TRUE, and a variable posMulti will be output; The variant site with multiple alleles can be analyzed using rareGWAMA.single.multiAllele function;}
}
\value{
A list of analysis results;
}
\description{
conditional analysis for single variant association test;
}
