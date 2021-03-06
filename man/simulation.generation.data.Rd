% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datagenertor.R
\name{simulation.generation.data}
\alias{simulation.generation.data}
\title{Generate simulation data}
\usage{
simulation.generation.data(n = n, coeff = coeff, matrix.option = 1,
  collinear = collinear, sig = sig, split.prop = 0.8, option = 1)
}
\arguments{
\item{n}{The number of rows in the \code{data.frame}}

\item{coeff}{A vector a true coefficients}

\item{matrix.option}{1: Using an Exchangeable correlation matrix to simulate the predictors\cr
2: Using an Autoregressive correlation matrix to simulate the predictors\cr}

\item{collinear}{The correlation levels within the \code{matrix.option}}

\item{sig}{The model inherent error, the \eqn{\sigma^2}}

\item{split.prop}{An element specifying the training proportion. Note the testing proportion will be 1 - the training proportion.}

\item{option}{1: split the dataset according to \code{c(split.prop, 1 - split.prop)}\cr
2: Use the whole dataset. Note When \code{option = 2}, the \code{split.prop} will be ignored\cr}
}
\value{
A list of elements:
         \item{ytrain}{return ONLY when \code{option = 1}, this is the training proportion for y}
         \item{xtrain}{return ONLY when \code{option = 1}, this is the training set for the predictors}
         \item{ytest}{return ONLY when \code{option = 1}, this is the testing proportion for y}
         \item{xtest}{return ONLY when \code{option = 1}, this is the testing proportion for the predictors}
         \item{resample.data}{return ONLY when \code{option = 2}, this is the whole dataset}
}
\description{
This function generates data for the simulation studies. This function is extremely adaptive and allows for many different data set ups.
The function offers variability on: dimension of the \code{data.frame}, user-defined coefficient vector, 2 types of correlation matrixs,
different correlation levels, varation of the model error sigma and the option to split the \code{data.frame} into a training and a testing set.
}
\section{Details}{
 This is the data generating function for the simulation functions: \code{\link{simulation.collinear}} and \code{\link{simulation.adlasso}}
}

\examples{
#100 observations, true coefficient vector is c(1,2,3,4), using the exchangeable correlation matrix,
#correlation level within the exchangeable matrix is 0.3, model sigma is 2, splitting
#the dataset into 0.8/0.2.
result <- simulation.generation.data(n = 100, coeff = c(1,2,3,4), matrix.option = 1,
collinear = 0.3, sig = 2,split.prop = 0.8, option = 1)

}
\author{
Mokyo Zhou
}
