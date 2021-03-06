% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resi_boot.R
\name{resi_boot}
\alias{resi_boot}
\title{Internal lassoenet function}
\usage{
resi_boot(x = x, y = y, B.rep = B.rep, significance = significance,
  best.lambda = best.lambda, best.coef = best.coef, bestgamma = bestgamma,
  ridgelambda = ridgelambda, parallel = parallel, method = 1,
  dataa = data, xx.indices = c(1, 2, 3))
}
\arguments{
\item{x}{A \code{model.matrix} of the predictors.\cr}

\item{y}{A vector of response values for the model fitting.\cr}

\item{B.rep}{The number of residual bootstrappings to do.}

\item{significance}{The significance level of the confidence intervals e.g. 100(1-\eqn{\alpha})\%.\cr}

\item{best.lambda}{The optimum lambda selected from either the "OLS" or "ridge" method.}

\item{best.coef}{A vector of coefficients from either the "OLS" or "ridge" method for constructing inital weights.\cr}

\item{bestgamma}{The optimum gamma selected from either the "OLS" or "ridge" method.\cr}

\item{ridgelambda}{The optimum lambda for the ridge regression (If the method "ridge" has been selecteed).\cr}

\item{parallel}{Parallelisation}

\item{method}{1 indicates "OLS", 2 indicates "ridge".}

\item{dataa}{Your full \code{data.frame}}

\item{xx.indices}{Locations of the predictors within dataa.}
}
\value{
The 100(1-\eqn{\alpha})\% confidence intervals for the point estimates from using residual bootstrapping.
}
\description{
Internal lassoenet function
}
\section{Details}{
 These are not intended for use by users. This function takes in outputs from one of \code{\link{OLS.approach}} or \code{\link{OLS.robust}} or \code{\link{ridge.approach}} or \code{\link{ridge.robust}} and
computes the corresponding 100(1-\eqn{\alpha})\% confidence intervals for the parameters. The return from this function will enter the funtion \code{\link{Adlasso}} for futher wrapping.
}

\author{
Mokyo Zhou
}
