% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/en.robust.split.R
\name{en.robust.split}
\alias{en.robust.split}
\title{Internal lassoenet function}
\usage{
en.robust.split(x = x, y = y, test = test, yv = yv,
  alpha.seq = alpha.seq, err.curves = err.curves,
  result.matrix = result.matrix, parallel = parallel)
}
\arguments{
\item{x}{A \code{model.matrix} of the preditors for the training data.\cr}

\item{y}{A vector of response values on the training set.\cr}

\item{test}{A \code{model.matrix} of the predictors for the testing data.\cr}

\item{yv}{A vector of response values on the testing set.\cr}

\item{alpha.seq}{The alpha sequence for the Elastic Net.\cr}

\item{err.curves}{The number of error curves to be fitted. Default is 0.\cr}

\item{result.matrix}{A matrix for storing results.\cr}

\item{parallel}{Parallelisation}
}
\value{
A vector of results for the best Elastic Net model along with some plots of the error curves, under the condition where the full dataset has been divided into a traning and a test set. The return from this function will enter \code{\link{prediction_ElasticNet}}.
}
\description{
Internal lassoenet function
}
\section{Details}{
 These are not intended for use by users. This function is one of the main engines for the Elastic Net computation. This function is used when the user would like to split the dataset into
 a training and a test set. This together with the function \code{\link{split_en}} forms the computation operator for the Elastic Net when using the a traning and a testing set. The return from this function will enter \code{\link{prediction_ElasticNet}} for futher wrapping.
}

\author{
Mokyo Zhou
}
