% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rareGWAMA.single.R
\name{re2.mr}
\alias{re2.mr}
\title{mixed effect meta-regression}
\usage{
re2.mr(b.vec, se.vec, af.pca)
}
\arguments{
\item{b.vec}{genetic effect estimates}

\item{se.vec}{standard deviation for the genetic effect estimates}

\item{af.pca}{the covariates to be included in the meta-regression}
}
\value{
a list consists of p.value, statistic, gamma.est and tau2;
}
\description{
mixed effect meta-regression
}
