% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rareGWAMA.conditional.single.R
\name{getCondUV}
\alias{getCondUV}
\title{get conditional U and V; implement our PCBS method for conditional association analyses;}
\usage{
getCondUV(...)
}
\arguments{
\item{ix.canddidate}{candidate variant indices;}

\item{ix.known}{conditioned variant indices}

\item{dat}{the list with contributed summary stat;}

\item{r2}{LD matrices;}

\item{sandwich}{use sandwich estimator or not; default is not to apply sandwich estimator;}
}
\value{
a list with condiitonal association analysis statistics;
}
\description{
get conditional U and V; implement our PCBS method for conditional association analyses;
}
