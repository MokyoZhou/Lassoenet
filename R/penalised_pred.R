#'Internal lassoenet functions
#'
#'Internal lassoenet functions
#'
#' @param data A well-cleaned \code{data.frame}.\cr
#' @param parallel parallelisation.\cr
#' @param response The location of the response within the \code{data.frame}.\cr
#' @param x.indices The locations of the predictors withint the \code{data.frame}.\cr
#' @param err.curves The number of error curves to fit.\cr
#' @param type.lambda Either "lambda.min" or "lambda.1se".
#' @return A vector of outputs and some plots related to the best model. The return from this function will enter the \code{\link{interactiver}} function.
#'
#' @section Details: These are not intended for use by users. This function is the overall wrapper for the prediction focused path. It takes returns from both the \code{\link{prediction_Lasso}} and \code{\link{prediction_ElasticNet}}
#' and put these through the comparison functions \code{\link{prediction.nonsplit.result}} and \code{\link{prediction.split.result}}. The return from this function will enter the \code{\link{interactiver}} function.
#'
#'
#' @author Mokyo Zhou
#'
#'@export
#'
#'


penalised_pred <- function(data = data, parallel = parallel, response = response, x.indices = x.indices, err.curves = err.curves, type.lambda = type.lambda){
  A10 <- readline("8. The algorithm will run both the Lasso and Elastic Net.
Please give the step size of the corresponding alpha sequence. Note the alpha
can only be ranged between [0,1], choose the step size appropreately. e.g. 0.1")
  Ans <- ifelse(is.na(A10),stop("Please enter a number for the alpha step size"),NUL <-0)
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(A10)),stop("please only input a numeric value for alpha step size"),NUL<-0))
  A10 <- as.numeric(A10)
  Ans <- ifelse(length(A10) != 1,stop("please make sure the step size for alpha is a single number"),NUL<-0)
  Ans <- ifelse(as.numeric(A10)<=0 || as.numeric(A10) >=1 ,stop("please make sure the step size is within the range (0,1)"),NUL<-0)


  cat("step size of alpha sequence =", A10)

  if (dim(data)[1]/(dim(data)[2]-1) < 10) {

    lasso.results.nosplit <- prediction_Lasso(data = data, interactive = TRUE, err.curves = err.curves, splits = c(0,0), x.indices = x.indices, response=response,type.lambda=type.lambda,parallel = parallel)


    elasticnet.results.nosplit <- prediction_ElasticNet(data = data, interactive = TRUE, err.curves = err.curves, splits = c(0,0),
                                                        x.indices = x.indices, response=response, step.size = A10, type.lambda = type.lambda,parallel = parallel)


    #comparsion here (non split)
    if(lasso.results.nosplit[3] < elasticnet.results.nosplit[4]){
      best.lasso <- prediction.nonsplit.result(best.lasso.result = lasso.results.nosplit, best.EN.result = elasticnet.results.nosplit,
                                 data = data, x.indices = x.indices, response=response, alph = 1, parallel = parallel)

      return(best.lasso)

    }else{
      best.ElasticNet <- prediction.nonsplit.result(best.lasso.result = lasso.results.nosplit, best.EN.result = elasticnet.results.nosplit,
                                data = data, x.indices = x.indices, response=response, alph = elasticnet.results.nosplit[2],parallel = parallel)

      return(best.ElasticNet)
    }


  }

  if (dim(data)[1]/(dim(data)[2]-1) > 10){

    cat("                                                                                                      ")

    A4 <- readline("9. Your dataset is relatively large, which means that you
can do a data split that could potentially give you a more relastic prediction
assessment. The drawback is that you will have less data to train with. Split Yes or No?")
    Ans <- ifelse(A4 != "Yes" & A4 != "No",stop("please only input /Yes/ or No/"),NUL<-0)
    cat("Splitting Yes or No =", A4)

    if (A4 == "Yes"){
      cat("                                                                                             ")

      A5 <- readline("10. please specify the splitting proportion with the first
number indicates the proportion of train and the second indicates the proportion
of test e.g. 0.7, 0.3")
      AA <- strsplit(A5, ",")
      Ans <- ifelse(length(AA[[1]])!=2,stop("please make sure a value for training and a value for the testing proportion are given"),NUL<-0)
      Ans <- suppressWarnings(ifelse(is.na(as.numeric(AA[[1]][1])),stop("please only input a numeric value for the training proportion"),NUL<-0))
      Ans <- suppressWarnings(ifelse(is.na(as.numeric(AA[[1]][2])),stop("please only input a numeric value for the testing proportion"),NUL<-0))
      splitpro <- c(as.numeric(AA[[1]][1]),as.numeric(AA[[1]][2]))
      Ans <- ifelse(sum(splitpro)!=1,stop("please make sure the proportions do sum to 1"),NUL<-0)
      Ans <- ifelse(sum(splitpro>=0)!=2,stop("please make sure both training and testing proportions are  > 0"),NUL<-0)
      Ans <- ifelse(sum(splitpro<=1)!=2,stop("please make sure the proportions are less than 1"),NUL<-0)
      cat(c("splitting proportion =", A5))

      cat("                                                                                              ")

      readline("we strongly advise you to try multiple different splitting
proportions and find one that gives the highest performance")

      lasso.results.split <- prediction_Lasso(data = data, interactive = TRUE, err.curves = err.curves, splits = splitpro, x.indices = x.indices, response=response,type.lambda=type.lambda,parallel = parallel)


      elasticnet.results.split <- prediction_ElasticNet(data = data, interactive =TRUE, err.curves = err.curves, splits = splitpro,
                                                        x.indices = x.indices, response=response, step.size = A10, type.lambda=type.lambda,parallel = parallel)

      cat("                                                                                               ")
      #comparsion here (split)
      A13 <- readline("11. please choose the preferred prediction metric for comparision:
MSE(Mean Squared Error) or MAE(Mean Absolute Error, please only enter either MSE or MAE")
      cat(c("prediction metric for comparison =", A13))
      Ans <- ifelse(A13 != "MSE" & A13 != "MAE",stop("please only input MSE or MAE"),NUL<-0)
      if (A13 == "MSE"){matric <- 3}; if(A13=="MAE"){matric <- 5}
      if(lasso.results.split[matric] < elasticnet.results.split[matric + 1]){
        best.lasso <- prediction.split.result(best.lasso.result = lasso.results.split, best.EN.result = elasticnet.results.split,
                                              data = data, x.indices = x.indices, response=response, alph = 1,
                                              splits = splitpro,parallel = parallel)

        return(best.lasso)

      }else{
        best.ElasticNet <- prediction.split.result(best.lasso.result = lasso.results.split, best.EN.result = elasticnet.results.split,
                                              data = data, x.indices = x.indices, response=response, alph = elasticnet.results.split[2],
                                              splits = splitpro,parallel = parallel)

        return(best.ElasticNet)
      }

    }else if (A4 == "No"){
      lasso.results.nosplit <- prediction_Lasso(data = data, interactive = TRUE, err.curves = err.curves, splits = c(0,0), x.indices = x.indices, response=response,type.lambda=type.lambda,parallel = parallel)


      elasticnet.results.nosplit <- prediction_ElasticNet(data = data, interactive = TRUE, err.curves = err.curves, splits = c(0,0),
                                                          x.indices = x.indices, response=response, step.size = A10, type.lambda = type.lambda,parallel = parallel)

      cat("                                                                                               ")
      cat("                                                                                               ")

      if(lasso.results.nosplit[3] < elasticnet.results.nosplit[4]){
        best.lasso <- prediction.nonsplit.result(best.lasso.result = lasso.results.nosplit, best.EN.result = elasticnet.results.nosplit,
                                                 data = data, x.indices = x.indices, response=response, alph = 1,parallel = parallel)

        return(best.lasso)

      }else{
        best.ElasticNet <- prediction.nonsplit.result(best.lasso.result = lasso.results.nosplit, best.EN.result = elasticnet.results.nosplit,
                                                      data = data, x.indices = x.indices, response=response, alph = elasticnet.results.nosplit[2],parallel = parallel)

        return(best.ElasticNet)
      }
    }

  }
}








