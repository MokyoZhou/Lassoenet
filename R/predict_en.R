#' Find the optimum pair (\eqn{\alpha}, \eqn{\lambda}) for an Elastic Net model
#'
#' This function takes in inputs defined by the user and computes the optimum \eqn{\alpha} and \eqn{\lambda} for an Elastic Net model. The function is very flexible and allows for many different settings
#' such as, data splitting, repeated error curves and a definable \eqn{\alpha} grid. This function also fully supports multiple-cores parallelisation. The main fitting process is cv.glmnet() from the
#' package glmnet.
#'
#' @param data A well-cleaned \code{data.frame} which will be used for modelling. The \code{data.frame} is also required to have more rows than columns.
#' @param x.indices The coordinates of the predictors that you would like to model with. Please provide a vecotr of locations e.g. seq(2,6).
#' @param response The location of the response within the \code{data.frame}.
#' @param err.curves Due to the fact that the cross-validation process is random, the result can vary qutie a bit (if without a seed). In order to stabilise the CV process, the function fits a collection of
#'                   Elastic Net models, each with a different \eqn{\alpha} value, multiple times. Therefore, for EACH \eqn{\alpha}, we will create multiple error curves over a range of \eqn{\lambda}s and the optimum
#'                   pair (\eqn{\alpha}, \eqn{\lambda}) is the pair that has the lowest averaged error curves value (local optimum). We have an optimum pair (\eqn{\alpha}, \eqn{\lambda}) for each \eqn{\alpha} and the global optimum pair
#'                   is the pair that has the overall lowest averaged error curved value.  Note, with this setting, the process tends to be slow. Thus, we highly suggest
#'                   multiple-cores parallelisation. Note you can set this argument to 0 if you do not wish to stabilise the process, in which case the seed (1234567) will be used for the CV process.
#'                   A postive integer indicates the stabilisation process is desired.
#'                   For more information about how this works, please see section details below.
#' @param splits A element specifying the training proportion and test proportion will be set as 1- traing.proportion. Note if you set \code{splits = 0},
#'               the function will use the whole dataset for modelling. The default is \code{splits = 0}.
#' @param step.size The step size of the \eqn{\alpha} grid. An adjustable value within [0,1]. Default is \code{step.size = 0.2}. Considering the time consumption, this should be chosen very carefully if \code{err.curves >0} is also desired.
#' @param type.lambda Either "lambda.min" or "lambda.1se". Default is "lambda.min. Note when \code{err.curves >0}, this argument will not be used.
#' @param interactive If you are running this function, please ALWAYS keep this argument to FALSE, which is the default.
#' @param parallel parallelisation supported,default is FALSE.
#' @return a list with elements:
#' \item{seed}{if \code{err.curves = 0}, the seed (1234567) will be used to compute \eqn{\alpha} and \eqn{\lambda}.\cr}
#' \item{number of err.curves}{if \code{err.curves >0}, this shows how many error curves there are (for each \eqn{\alpha}).\cr}
#' \item{best alpha}{if \code{err.curves = 0} (no stabilization), this is a part of the optimum pair (\eqn{\alpha}, \eqn{\lambda}) that has the lowest cross validation error (out of a single 2D grid search) with seed (1234567). This is the usual way of finding out (\eqn{\alpha}, \eqn{\lambda}).\cr
#'                   if \code{err.curves} > 0 (with stabilization), this is a part of the global optimum pair (\eqn{\alpha}, \eqn{\lambda}) that has the overall lowest averaged error curves value among all the local optimum pairs.\cr}
#' \item{best lambda}{A part of the global optimum pair (\eqn{\alpha}, \eqn{\lambda}). For more information on how the \eqn{\alpha} and \eqn{\lambda} get selected, please see the detail section.\cr}
#' \item{prediction error}{if \code{err.curves = 0} and \code{splits = 0}, (no error curves and no splitting), this is the cross validation score associated with the best (\eqn{\alpha},\eqn{\lambda}) pair from the a 2D grid search with seed (1234567).\cr
#' if \code{err.curves > 0} and \code{splits = 0}, (error curves but no splitting), this is the overall lowest averaged cross-validation scores associated with the global optimum (\eqn{\alpha}, \eqn{\lambda}) pair from using the whole dataset.\cr
#' if \code{err.curves = 0} and \code{splits != 0}, (no error curves but with split), this is the out of sample Mean-squared error on the test set with the optimum pair (\eqn{\alpha},\eqn{\lambda}) selected with a 2D grid search with seed (1234567) on the training set.\cr
#' if \code{err.curves > 0} and \code{splits != 0}, (error curves and splitting), this is the out of sample Mean-squared error on the test set with the global optimum pair (\eqn{\alpha},\eqn{\lambda}) selected by averaging across the error curves on the training set.\cr}
#' \item{prediction_lower}{The function will only compute this when \code{splits = 0}. For more information see details below.\cr}
#' \item{prediction_upper}{The function will only compute this when \code{splits = 0}, For more information see details below.\cr}
#' \item{rooted prediction error}{The function will only compute this when \code{splits != 0}. This is the out of sample root mean-square error on the test set.\cr}
#' \item{absolute prediction error}{The function will only compute this when \code{splits != 0}. This is the out of sample mean absolute error on the test set.\cr}
#'
#' @section Details: This function further develops on the cv.glmnet() function from the glmnet package to allow for more flexibility.
#' The function allows users the option to split the dataset into a traning set and a test set, which usually
#' gives more realistic assessment of perdictive performance than using the whole dataset.
#'
#' The function also offers an alternative to compute the global optimum (\eqn{\alpha}, \eqn{\lambda}) pair by averaging across the error curves instead of using a fixed seed. More specifically, this stabilisation
#' takes a double looping structure where the outer layer contains the \eqn{\alpha} grid and the inner loop builds multiple Elastic Net models for each \eqn{\alpha} in the outer layer.
#' In this way, for each \eqn{\alpha}, we create say, B, Elastic Net models and thus, this results in B error curves over a range of \eqn{\lambda}s. For each \eqn{\alpha} then, the
#' function finds the local optimum pair (\eqn{\alpha}, \eqn{\lambda}) by averaging across these error curves and find the pair that has the lowest averaged cross validation errors. After the function finds
#' all the local optimum pairs, the golbal optimum pair is the pair that has the overall lowest averaged cross validation error. From experneice, for medium
#' size datasets, with \code{err.curves} larger than 1500, the global optimum (\eqn{\alpha}, \eqn{\lambda}) will usually converge to stable values that consistently
#' achieves the overall lowest averaged across error curves value.
#'
#' When \code{splits = 0} (no spliting), the 95 percent confidence interval for such a prediction error (overall lowest averaged error curves value) is generated as follows:
#' from the corresponding error curves for that optimum \eqn{\alpha}, the cross-validation scores corresponding to the global optimum (\eqn{\alpha}, \eqn{\lambda}) are extracted and the command quantile()
#' is then used to compute the 95 percent confidence interval. Do not worry, we will also provide you with a plot that contains these error curves from the output, so you can
#' see how it is that the optimum pair values got selected and with its 95 percent CI around it. When we are not stabilising the process e.g.
#' \code{err.curve = 0} but still using the whole dataset \code{splits = 0}, we compute the (\eqn{\alpha}, \eqn{\lambda}) pair with the seed (1234567) and the associated CI is computed by using the standard error
#' provided by the glmnet package and assuming normality. When we do split the dataset however, we provide out of sample mean squared error(MSE), RMSE and MAE instead of 95 percent CI.
#'
#'
#'
#' @author Mokyo Zhou
#'
#' @examples
#'
#' library(glmnet)
#' data(QuickStartExample)
#' #please NOTE: you can access "QuickStartExample" by using data.frame(y,x)
#'
#'
#' #non-split, no error curves, step.size =0.2, using lambda.min
#' result <- prediction_ElasticNet(data = data.frame(y,x), x.indices = seq(2,21),
#' response = 1, err.curves = 0 ,splits = 0, step.size =0.2,
#' type.lambda="lambda.min")
#'
#' #0.8 /0.2 split, no error curves, step.size=0.1, using lambda.1se
#' result <- prediction_ElasticNet(data = data.frame(y,x), x.indices = seq(2,21),
#' response = 1, err.curves = 0, splits = 0.8, step.size = 0.1,
#' type.lambda = "lambda.1se")
#'
#' #non-split, but with 100 error curves, step.size=0.3, with parallel (2 cores)
#' #cl <- parallel::makeCluster(2)
#' #doParallel::registerDoParallel(cl)
#' result <- prediction_ElasticNet(data = data.frame(y,x), x.indices = seq(2,21),
#' response = 1, err.curves = 100, splits = 0, step.size = 0.3, parallel = TRUE)
#'
#' #0.8 / 0.2 split, with 100 error curves, step.size=0.2, with parallel (2 cores)
#' #cl <- parallel::makeCluster(2)
#' #doParallel::registerDoParallel(cl)
#' result <- prediction_ElasticNet(data = data.frame(y,x), x.indices = seq(2,21),
#' response = 1, err.curves = 100, splits = 0.8, step.size = 0.2, parallel = TRUE)
#'
#'
#' @export
#'
#'
#'
#'



prediction_ElasticNet <- function(data = data, x.indices = x.indices, response = response, err.curves = 0 ,
                                  splits = 0, step.size =0.2, type.lambda="lambda.min", interactive = FALSE, parallel = FALSE){



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
  Ans <- ifelse(is.numeric(err.curves)==FALSE,stop("Please make sure you have a numeric input for err.curves"),NUL<-0)
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(err.curves)),stop("please only input a single value for the number of error curves"),NUL<-0))
  Ans <- ifelse(err.curves != floor(err.curves) ,stop("Please make sure the number of error curves is a whole number"),NUL<-0)
  Ans <- ifelse(err.curves < 0 ,stop("Please make sure the number of error curves is positive"),NUL<-0)
  Ans <- ifelse(length(err.curves)!=1 ,stop("Please make sure only enter a single number for the number of error curves"),NUL<-0)
  if(interactive==FALSE){
    Ans <- ifelse(is.na(as.numeric(splits)),stop("please only enter a number for the training proportion"),NUL<-0)
    Ans <- ifelse(length(splits)!=1,stop("please only enter a number for the training proportion"),NUL<-0)
    Ans <- ifelse(splits < 0 || splits >1, stop("please make sure the training proportion is within [0,1]"),NUL<-0)
  }
  Ans <- ifelse(is.na(step.size),stop("please make sure the step size for alpha is a number"),NUL<-0)
  Ans <- ifelse(length(step.size) != 1,stop("please make sure the step size for alpha is a single number"),NUL<-0)
  Ans <- ifelse(step.size <= 0 || step.size >=1,stop("please make sure the alpha step size is within [0,1]"),NUL<-0)
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
    no.split.en <- no_split_en(err.curves = err.curves , x = data[[2]], y = data[[1]], step.size = step.size, type.lambda=type.lambda, parallel = parallel)
    if(interactive == FALSE & err.curves ==0){
      no.split.en <- list("seed" = no.split.en[1],"best alpha" = no.split.en[2], "best lambda" = no.split.en[3],
                          "prediction error" = no.split.en[4],"prediciton_lower" = no.split.en[5], "prediction_upper" = no.split.en[6])
    }
    if(interactive ==FALSE & err.curves > 0){
      no.split.en <- list("number of err.curves" = no.split.en[1],"best alpha" = no.split.en[2], "best lambda" = no.split.en[3],
                             "prediction error" = no.split.en[4], "prediciton_lower" = no.split.en[5], "prediction_upper" = no.split.en[6])
    }
    return(no.split.en)


  }else if(sum(splits)!=0){
    splitdata <- datasplit(data = data, splits = splits)
    train.data <- convertor(splitdata[[1]], response = response, x.indices = x.indices)
    testt <- convertor(splitdata[[2]], response = response, x.indices = x.indices)
    split.en <- split_en(err.curves = err.curves, x.train = train.data[[2]], y.train = train.data[[1]], x.test = testt[[2]]
                         ,y.test = testt[[1]],step.size=step.size, type.lambda=type.lambda, parallel = parallel)
    if(interactive == FALSE & err.curves ==0){
      split.en <- list("seed" = split.en[1],"best alpha" = split.en[2], "best lambda" = split.en[3],
                          "prediction error" = split.en[4],"rooted prediciton error" = split.en[5], "absolute prediction error" = split.en[6])
    }
    if(interactive ==FALSE & err.curves > 0){
      split.en <- list("number of err.curves" = split.en[1],"best alpha" = split.en[2], "best lambda" = split.en[3],
                             "prediction error" = split.en[4], "rooted prediciton error" = split.en[5], "absolute prediction error" = split.en[6])
    }
    return(split.en)
  }




}



