#' Find the optimum \eqn{\lambda} for a Lasso model
#'
#' This function takes in inputs defined by the user and computes the optimum \eqn{\lambda} for a Lasso model. The function is very flexible and allows for many different settings
#' such as, data splitting and repeated error curves. This function also fully supports multiple-cores parallelisation. The main fitting process is cv.glmnet() from the
#'  package glmnet.
#'
#' @param data A well-cleaned \code{data.frame} which will be used for modelling. The \code{data.frame} is also required to have more rows than columns.
#' @param x.indices The coordinates of the predictors that you would like to model with. Please provide a vecotr of locations e.g. seq(2,6).
#' @param response The location of the response within the \code{data.frame}.
#' @param err.curves Due to the fact that the Cross-validation process is random, it is very likely that the result will vary quite a bit (if without a seed).
#'                   This function offers to fit the model multiple times (thus, creating multiple error curves over a range of \eqn{\lambda}s) and average across
#'                   these multiple error curves. The optimum \eqn{\lambda} within the range is the one that has the lowest averaged error. Note you can set this argument
#'                   to 0 if you do not wish to stabilise the process, in which case the seed (1234567) will be used for the CV process. A positive integer indicates
#'                   the number of error curves to be fitted. Default is 0.
#' @param splits A element specifying the training proportion and the test proportion will be set as 1 - traning.proportion. Note if you set \code{splits = 0},
#'               the function will use the whole dataset for modelling. The default is \code{splits = 0}.
#' @param type.lambda Either "lambda.min" or "lambda.1se". Default is "lambda.min. Note when \code{err.curves >0}, this argument will not be used.
#' @param interactive If you are running this function, please ALWAYS keep this argument to FALSE, which is the default.
#' @param parallel parallelisation supported,default is FALSE.
#' @return a list with elements:
#' \item{seed}{if \code{err.curves = 0}, the seed (1234567) will be used to compute the \eqn{\lambda}.\cr}
#' \item{number of err.curves}{if \code{err.curves > 0}, this argument shows how many error curves there are.\cr}
#' \item{best lambda}{if \code{err.curves = 0}, this is the \eqn{\lambda} that has the lowest cross validation error from a Lasso fit with seed (1234567).This is the usual way to select \eqn{\lambda}.\cr
#'  if \code{err.curves} > 0, this is the the \eqn{\lambda} that has the lowest averaged error curves value.\cr}
#' \item{prediction error}{if \code{err.curves = 0} and \code{splits = 0}, (no error curves and no splitting), this simply is the cross validation score (associated with the best \eqn{\lambda}) from the seed (1234567).\cr
#' if \code{err.curves > 0} and \code{splits = 0}, (error curves but no splitting), this is the averaged cross-validation scores associated with the optimum \eqn{\lambda} from using the whole dataset.\cr
#' if \code{err.curves = 0} and \code{splits != 0}, (no error curves but with split), this is the out of sample Mean-squared error on the test set with \eqn{\lambda} selected with seed(1234567) on the training set.\cr
#' if \code{err.curves > 0} and \code{splits != 0}, (error curves and splitting), this is the out of sample Mean-squared error on the test set with \eqn{\lambda} selected by averaging across the error curves on the training set.\cr}
#' \item{prediction_lower}{The function will only compute this when \code{splits = 0}. For more information see details below.\cr}
#' \item{prediction_upper}{The function will only compute this when \code{splits = 0}, For more information see details below.\cr}
#' \item{rooted prediction error}{The function will only compute this when \code{splits != 0}. This is the out of sample root mean-square error on the test set.\cr}
#' \item{absolute prediction error}{The function will only compute this when \code{splits != 0}. This is the out of sample mean absolute error on the test set.\cr}
#'
#' @section Details: This function further develops on the cv.glmnet() function from the glmnet package to allow for more flexibility.
#' More specifically, it allows users the option to split the dataset into a traning set and a test set, which usually
#' gives give more realistic assessment of perdictive performance than using the whole dataset.\cr
#'
#' The function also offers an alternative
#' to compute the optimum \eqn{\lambda} by averaging across the error curves instead of using a fixed seed. From experneice, for medium
#' size datasets, with \code{err.curves} larger than 1000, the optimum \eqn{\lambda} will usually converge to a stable value that consistently
#' achieves the lowest averaged across error curves value.\cr
#'
#' The 95 percent confidence interval for this lowest averaged error curves value,
#' when using the whole dataset e.g. \code{splits = 0}, is computed by using the quantile() command on the cross-validation scores
#' associated with this optimum \eqn{\lambda}. Do not worry, we will also provide you with a plot that contains these error curves from the output, so you can
#' see how it is that the optimum \eqn{\lambda} value got selected and with its 95 percent CI around it. When we are not stabilising the process e.g.
#' \code{err.curve = 0} but still using the whole dataset \code{splits = 0}, we compute the \eqn{\lambda} with the seed (1234567) and the associated CI is computed by using the standard error
#' provided by the glmnet package and assuming normality. When we do split the dataset however, we provide out of sample mean squared error(MSE), RMSE and MAE instead of 95 percent CI.
#'
#' @author Mokyo Zhou
#'
#' @examples
#'
#' library(glmnet)
#' data(QuickStartExample)
#' #please NOTE: you can access "QuickStartExample" by using: data.frame(y,x)
#'
#'
#' #non-split, no error curves, using lambda.min
#' result <- prediction_Lasso(data = data.frame(y,x), x.indices = seq(2,21),
#' response = 1, err.curves = 0,splits = 0)
#'
#' #0.8 /0.2 split, no error curves, using lambda.1se
#' result <- prediction_Lasso(data = data.frame(y,x), x.indices = seq(2,21),
#' response = 1, err.curves = 0, splits = 0.8,type.lambda = "lambda.1se")
#'
#' #non-split, but with 100 error curves with parallel (2 cores)
#' #cl <- parallel::makeCluster(2)
#' #doParallel::registerDoParallel(cl)
#' result <- prediction_Lasso(data = data.frame(y,x), x.indices = seq(2,21),
#' response = 1, err.curves = 100, splits = 0,parallel = TRUE)
#'
#' #0.8 / 0.2 split, with 100 error curves with parallel (2 cores)
#' #cl <- parallel::makeCluster(2)
#' #doParallel::registerDoParallel(cl)
#' result <- prediction_Lasso(data = data.frame(y,x), x.indices = seq(2,21),
#' response = 1, err.curves = 100, splits = 0.8,parallel = TRUE)
#'
#' @export
#'
#'
#'
#'


prediction_Lasso <- function(data = data, x.indices = x.indices, response = response, err.curves = 0,
                             splits = 0, type.lambda="lambda.min",interactive = FALSE, parallel = FALSE){



  if(dim(data)[2] > dim(data)[1]){
    stop("please make sure your dataset has more rows than columns")
  }
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(response)),stop("please only input a single value for the response coordinate"),NUL<-0))
  Ans <- ifelse((as.numeric(response) %in% seq(1,ncol(data)))==FALSE ,stop("please make sure you have the correct y coordinate"),NUL<-0)
  Ans <- ifelse(length(response)!=1 ,stop("please make sure you have only entered a single location for y"),NUL<-0)
  Ans <- ifelse(is.numeric(x.indices)==FALSE,stop("Please make sure you have numeric inputs for x.indices"),NUL<-0)
  Ans <- ifelse(sum((as.numeric(x.indices) %in% seq(1,ncol(data)))==FALSE)!=0 ,stop("please make sure you have the correct x coordinates"),NUL<-0)
  Ans <- ifelse(sum(as.numeric(x.indices) == as.numeric(response))!=0 ,stop("Please make sure Xs do not have the same coordinate as y"),NUL<-0)
  Ans <- ifelse(length(x.indices)!= length(unique(x.indices)) ,stop("Please make sure you have only entered each X coordinate once"),NUL<-0)
  Ans <- ifelse(length(x.indices)== 1,stop("Please make sure you have at least two X coordinates"),NUL<-0)
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(err.curves)),stop("please only input a single value for the number of error curves"),NUL<-0))
  Ans <- ifelse(is.numeric(err.curves)==FALSE,stop("Please make sure you have a numeric input for err.curves"),NUL<-0)
  Ans <- ifelse(err.curves != floor(err.curves) ,stop("Please make sure the number of error curves is a whole number"),NUL<-0)
  Ans <- ifelse((err.curves < 0) ,stop("Please make sure the number of error curves is positive"),NUL<-0)
  Ans <- ifelse(length(err.curves)!=1 ,stop("Please make sure only enter a single number for the number of error curves"),NUL<-0)
  Ans <- ifelse(is.na(as.numeric(splits)),stop("please only enter a number for the training proportion"),NUL<-0)
  if(interactive == FALSE){
    Ans <- ifelse(length(splits)!=1,stop("please only enter a number for the training proportion"),NUL<-0)
    Ans <- ifelse(splits < 0 || splits >1,stop("please make sure the training proportion is within [0,1]"),NUL<-0)
  }

  Ans <- ifelse(type.lambda != "lambda.min" & type.lambda != "lambda.1se" ,stop("Please make sure you have entered the correct lambda type"),NUL<-0)


  if(interactive == FALSE){
    if(splits != 0){
      splits <- c(splits, 1- splits)

    }else{
      splits <- c(0,0)
    }
  }

  if(sum(splits)==0){
    data <- convertor(data = data, response = response, x.indices = x.indices)
    no.split.lasso <- no_split_lasso(err.curves = err.curves , x = data[[2]], y = data[[1]], type.lambda=type.lambda,parallel = parallel)
    if(interactive == FALSE & err.curves ==0){
      no.split.lasso <- list("seed" = no.split.lasso[1],"best lambda" = no.split.lasso[2], "prediction error" = no.split.lasso[3],
                             "prediciton_lower" = no.split.lasso[4], "prediction_upper" = no.split.lasso[5])
    }
    if(interactive ==FALSE & err.curves > 0){
      no.split.lasso <- list("number of err.curves" = no.split.lasso[1],"best lambda" = no.split.lasso[2], "prediction error" = no.split.lasso[3],
                             "prediciton_lower" = no.split.lasso[4], "prediction_upper" = no.split.lasso[5])
    }

    return(no.split.lasso)


  }else if(sum(splits)!=0){
    splitdata <- datasplit(data = data, splits = splits)
    train.data <- convertor(splitdata[[1]], response = response, x.indices = x.indices)
    testt <- convertor(splitdata[[2]], response = response, x.indices = x.indices)
    split.lasso <- split_lasso(x.train = train.data[[2]], y.train = train.data[[1]], x.test = testt[[2]], y.test = testt[[1]],
                               err.curves = err.curves,type.lambda=type.lambda,parallel = parallel)

    if(interactive == FALSE & err.curves ==0){
      split.lasso <- list("seed" = split.lasso[1],"best lambda" = split.lasso[2], "prediction error" = split.lasso[3],
                             "rooted prediciton error" = split.lasso[4], "absolute prediction error" = split.lasso[5])
    }
    if(interactive ==FALSE & err.curves > 0){
      split.lasso <- list("number of err.curves" = split.lasso[1],"best lambda" = split.lasso[2], "prediction error" = split.lasso[3],
                          "rooted prediciton error" = split.lasso[4], "absolute prediction error" = split.lasso[5])
    }


    return(split.lasso)
  }




}


