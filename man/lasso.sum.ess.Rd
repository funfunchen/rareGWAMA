% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rareGWAMA.pred.R, R/rareGWAMA.twas.pred.R
\name{lasso.sum.ess}
\alias{lasso.sum.ess}
\title{lassosum with essential predictors;}
\usage{
lasso.sum.ess(b.vec, s.vec, r2.mat, n, group, lambda, alpha)

lasso.sum.ess(b.vec, s.vec, r2.mat, n, group, lambda, alpha)
}
\arguments{
\item{b.vec}{marginal effect of each genetic variant; note this has to be marginal effects;}

\item{s.vec}{the standard error for each genetic variant;}

\item{r2.mat}{the LD matrix for the predictors, which can be approximated from a reference panel;}

\item{n}{sample size;}

\item{group}{a group of categories for the predictors, should be of values 1,2, ... etc;}

\item{lambda}{tuning parameters for L2 penalties}

\item{alpha}{tuning parameters for L1 panelties}
}
\description{
lassosum with essential predictors;

lasso with summary statisitcs and multi penalties
}
