#' Generate simulation data
#'
#' This function generates data for the simulation studies. This function is extremely adaptive and allows for many different data set ups.
#' The function offers variability on: dimension of the \code{data.frame}, user-defined coefficient vector, 2 types of correlation matrixs,
#' different correlation levels, varation of the model error sigma and the option to split the \code{data.frame} into a training and a testing set.
#'
#' @param n The number of rows in the \code{data.frame}
#' @param coeff A vector a true coefficients
#' @param matrix.option 1: Using an Exchangeable correlation matrix to simulate the predictors\cr
#'                      2: Using an Autoregressive correlation matrix to simulate the predictors\cr
#' @param collinear The correlation levels within the \code{matrix.option}
#' @param sig The model inherent error, the \eqn{\sigma^2}
#' @param split.prop An element specifying the training proportion. Note the testing proportion will be 1 - the training proportion.
#' @param option 1: split the dataset according to \code{c(split.prop, 1 - split.prop)}\cr
#'               2: Use the whole dataset. Note When \code{option = 2}, the \code{split.prop} will be ignored\cr
#' @return A list of elements:
#'          \item{ytrain}{return ONLY when \code{option = 1}, this is the training proportion for y}
#'          \item{xtrain}{return ONLY when \code{option = 1}, this is the training set for the predictors}
#'          \item{ytest}{return ONLY when \code{option = 1}, this is the testing proportion for y}
#'          \item{xtest}{return ONLY when \code{option = 1}, this is the testing proportion for the predictors}
#'          \item{resample.data}{return ONLY when \code{option = 2}, this is the whole dataset}
#'
#' @section Details: This is the data generating function for the simulation functions: \code{\link{simulation.collinear}} and \code{\link{simulation.adlasso}}
#'
#' @author Mokyo Zhou
#'
#' @examples
#'#100 observations, true coefficient vector is c(1,2,3,4), using the exchangeable correlation matrix,
#'#correlation level within the exchangeable matrix is 0.3, model sigma is 2, splitting
#'#the dataset into 0.8/0.2.
#' result <- simulation.generation.data(n = 100, coeff = c(1,2,3,4), matrix.option = 1,
#' collinear = 0.3, sig = 2,split.prop = 0.8, option = 1)
#'
#' @export
#'



simulation.generation.data<- function(n = n, coeff = coeff, matrix.option = 1, collinear = collinear, sig = sig,
                                      split.prop = 0.8, option = 1){


  Ans <- ifelse(is.na(as.numeric(n)),stop("please make sure the number of observation is a number"),NUL<-0)
  Ans <- ifelse(n < 0,stop("please make sure the number of observation is a positive"),NUL<-0)
  Ans <- ifelse(n != floor(n),stop("please make sure the number of observation is a whole number"),NUL<-0)
  Ans <- ifelse(dim(n)[1] < 5 * length(coeff),stop("please make sure you dataset has observations at least 5 times the length of your coeff"),NUL<-0)
  Ans <- ifelse(is.na(as.numeric(n)),stop("please make sure the number of observation is a number"),NUL<-0)
  Ans <- ifelse(length(coeff) ==1,stop("Please make sure you have more than 1 element within the coeff"),NUL<-0)
  Ans <- ifelse(sum(is.na(as.numeric(coeff)))!=0,stop("Please make sure all the elements within coeff are numbers"),NUL<-0)
  Ans <- ifelse(is.na(as.numeric(collinear)),stop("please make sure the collinear is a number"),NUL<-0)
  Ans <- ifelse(length(collinear) != 1,stop("Please make sure you have only entered a single number for collienar"),NUL<-0)
  Ans <- ifelse(collinear < -1 || collinear > 1,stop("Please make sure collinear is within [-1,1]"),NUL<-0)
  Ans <- ifelse(is.na(as.numeric(sig)),stop("please make sure the sig is a number"),NUL<-0)
  Ans <- ifelse(length(sig) != 1,stop("Please make sure you have only entered a single number for sig"),NUL<-0)
  Ans <- ifelse(sig < 0 ,stop("Please make sure the sig is a positive number"),NUL<-0)
  Ans <- ifelse(is.na(as.numeric(split.prop)),stop("please only enter a number for the training proportion"),NUL<-0)
  Ans <- ifelse(length(split.prop)!=1,stop("please only enter a number for the training proportion"),NUL<-0)
  Ans <- ifelse(split.prop < 0 || split.prop >1, stop("please make sure the training proportion is within [0,1]"),NUL<-0)
  Ans <- ifelse(option != 1 & option != 2,stop("Please make sure you have the correct option input, 1:split, 2: nonsplit"),NUL<-0)
  Ans <- ifelse(matrix.option !=1 & matrix.option != 2,stop("Please make sure you have the correct matrix.option input, 1:exchangeable, 2:autoregressive"),NUL<-0)

  p <- length(coeff)

  #Defining levels of collinearity between the columns
  if(matrix.option != 1){
    CovMatrix <- outer(1:p, 1:p, function(x,y) {collinear^abs(x-y)})
  }else{
    CovMatrix <- matrix(collinear,p,p)
    diag(CovMatrix) <- 1
  }

  #undo the inherent seed 1234567
  set.seed(Sys.time())

  #Obtaining X data
  X <- MASS::mvrnorm(n, rep(0,p), CovMatrix)

  #Obtaining Y
  Y <- X %*% as.matrix(coeff) + stats::rnorm(n,0,sig)

  #resampled data
  resample.data <- data.frame(Y, X)

  #data splitting
  if(option == 1){
    TrainTest <- datasplit(resample.data, splits = c(split.prop, 1 - split.prop))

    #Train
    Train <- TrainTest[[1]]; xtrain <- as.matrix(Train[,2:ncol(Train)]); ytrain <- Train[,1]

    #Test
    Test <- TrainTest[[2]]; xtest <- as.matrix(Test[,2:ncol(Test)]); ytest <- Test[,1]


    return(list("ytrain"=ytrain,"xtrain" = xtrain,"ytest"= ytest,"xtest"=xtest))

  }

  #whole data return
  return(resample.data)

}



