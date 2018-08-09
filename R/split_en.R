#' Internal lassoenet function
#'
#' Internal lassoenet function
#'
#' @param err.curves The number of error curves to be fitted. Default is 0.\cr
#' @param x.train A \code{model.matrix} of the preditors for the training data.\cr
#' @param y.train A vector of response values on the training set.\cr
#' @param x.test A \code{model.matrix} of the predictors for the testing data.\cr
#' @param y.test A vector of response values on the testing set.\cr
#' @param step.size The step size of the alpha sequence.\cr
#' @param type.lambda Either "lambda.min" or "lambda.1se", default is "lambda.min.\cr
#' @param parallel Parallelisation
#' @return A vector of results for the best Elastic Net model under the condition where the full dataset has been splitted into a training and a testing set. The return from this function will enter \code{\link{prediction_ElasticNet}}.
#'
#' @section Details: These are not intended for use by users. This function is one of the main engines for the Elastic Net computation. This function is used when the user would like to split the dataset into
#'  a training and a test set. This together with the function \code{\link{en.robust.split}} forms the computation operator for the Elastic Net when using the a traning and a testing set. The return from this function will enter \code{\link{prediction_ElasticNet}} for futher wrapping.
#'
#'
#'
#' @author Mokyo Zhou
#'
#'
#'
#' @export
#'
#'
#'




split_en<- function (err.curves = 0, x.train = x.train, y.train=y.train, x.test = x.test, y.test = y.test,
                     step.size=0, type.lambda=type.lambda,parallel = parallel){
  alpha.seq <-seq(0,1,step.size)
  ifelse(alpha.seq[length(alpha.seq)]==1, alpha.seq <- alpha.seq[c(-1,-length(alpha.seq))],
         alpha.seq <- alpha.seq[c(-1)])


  i <- NULL
  if(err.curves==0){
    type.lambda <- type.lambda


    z <- ifelse(type.lambda=="lambda.min",z <- 9, z <- 10)
    `%dopar%` <- foreach::`%dopar%`


    result.matrix <-foreach::foreach(i = alpha.seq, .combine = rbind) %dopar% {

      set.seed(1234567);cv <-glmnet::cv.glmnet(x.train,y.train, parallel = parallel, alpha = i)
      prediction_error <- cv$cvm[which(cv$lambda==cv[z])]
      data.frame(i, cv[z], prediction_error)
    }
    result.matrix <- as.matrix(result.matrix)


    optimum_pair <-result.matrix[which(result.matrix[,3]== min(result.matrix[,3])),]

    best.en <- glmnet::glmnet(x.train,y.train,alpha=optimum_pair[1],lambda=optimum_pair[2])

    predictions <- predict(best.en, newx = x.test)

    #MSE
    prediction_error <- sum((y.test - predictions)^2)/nrow(x.test)
    #Root MSE
    root_prediction_error <- sqrt(prediction_error)
    #MAE
    absolute_error <- sum((abs(y.test-predictions)))/nrow(x.test)

    #refitting for number of nonzeros
    nonzero <- length(which(as.matrix(coef(best.en,s=optimum_pair[2]))[,1]!=0))-1

    split.en.vector <-c(1234567, optimum_pair[1],optimum_pair[2],prediction_error,root_prediction_error,absolute_error,nonzero)
    split.en.vector <- base::unname(split.en.vector)
    return(split.en.vector)

  }else if(err.curves!=0){
    result.matrix <-matrix(0,length(alpha.seq),3)
    best.en.result <- en.robust.split(x = x.train, y = y.train, test = x.test, yv = y.test, alpha.seq = alpha.seq, err.curves=err.curves,
                                      result.matrix = result.matrix, parallel = parallel)
    return(best.en.result)

  }


}

