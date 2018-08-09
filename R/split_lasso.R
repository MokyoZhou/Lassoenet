#' Internal lassoenet function
#'
#' Internal lassoenet function
#'
#' @param x.train A \code{model.matrix} of the preditors for the training data.\cr
#' @param y.train A vector of response values on the training set.\cr
#' @param x.test A \code{model.matrix} of the predictors for the testing data.\cr
#' @param y.test A vector of response values on the testing set.\cr
#' @param err.curves The number of error curves to be fitted. Default is 0.\cr
#' @param type.lambda Either "lambda.min" or "lambda.1se", default is "lambda.min.\cr
#' @param parallel Parallelisation
#' @return A vector of results for the best Lasso model under the condition where the training set has been used to train and tune, and the testing set is used for obtaining the out of sample prediction error. The return from this function will enter \code{\link{prediction_Lasso}}.
#'
#' @section Details: These are not intended for use by users. This function is one of the main engines for the Lasso computation. This function is used when the user would like to split the dataset into
#'  a training and a test set. This together with the function \code{\link{no_split_lasso}} form the overall computation operator for the Lasso models. The return from this function will enter \code{\link{prediction_Lasso}} for futher wrapping.
#'
#'
#'
#' @author Mokyo Zhou
#'
#'
#'
#' @export
#'








split_lasso <- function(x.train = x.train, y.train = y.train, x.test = x.test, y.test = y.test,
                        err.curves = 0,type.lambda=type.lambda,parallel = parallel){


  if(err.curves==0){

    set.seed(1234567);best.lambda<-glmnet::cv.glmnet(x.train, y.train,parallel = parallel)

    new.data <- x.test
    y.validate <- y.test
    y.validate <- as.numeric(y.validate)

    ifelse(type.lambda=="lambda.min",bestLambda <- best.lambda$lambda.min, bestLambda <- best.lambda$lambda.1se)
    best.lasso <- glmnet::glmnet(x.train,y.train,lambda = bestLambda)

    predictions <- predict(best.lasso, newx = new.data)


    #MSE
    prediction_error <- sum((y.validate - predictions)^2)/nrow(new.data)
    #Root MSE
    root_prediction_error <- sqrt(prediction_error)
    #MAE
    absolute_error <- sum((abs(y.validate-predictions)))/nrow(new.data)


    nonzero <- length(which(as.matrix(coef(best.lambda,s=bestLambda))[,1]!=0))-1

    split.lasso.vector <- c(1234567,bestLambda,prediction_error,root_prediction_error,absolute_error, nonzero)
    return(split.lasso.vector)


  }else{
    cat("                                                                                               ")
    #breaking the inherent seed (12345678) from `datasplit()`
    set.seed(Sys.time())
    `%dopar%` <- foreach::`%dopar%`
    lambdas <- foreach::foreach(i=1:err.curves, .combine = rbind) %dopar%{
      fit <- glmnet::cv.glmnet(x.train,y.train,parallel=parallel)
      data.frame(fit$lambda, fit$cvm)
    }


    lambda.full <- lambdas
    lambdas <- aggregate(lambdas[, 2], list(lambdas$fit.lambda), mean)

    plot(log(lambda.full[,1]),lambda.full[,2],pch=20,col="black",main="Averaging across the error curves (Lasso, Training set)",
         xlab = "log lambdas",ylab="MSE")
    lines(log(lambdas[,1]),lambdas[,2],col="blue")

    # select the best one
    bestindex <- which(lambdas[2]==min(lambdas[2]))
    bestlambda <- lambdas[bestindex,1]
    abline(v=log(bestlambda),lty=2)


    lasso <- glmnet::glmnet(x.train,y.train)

    new.data <- x.test
    y.validate <- y.test
    y.validate <- as.numeric((y.validate))

    predictions <- predict(lasso, newx = new.data, s = bestlambda)

    #MSE
    prediction_error <- sum((y.validate - predictions)^2)/nrow(new.data)
    #Root MSE
    root_prediction_error <- sqrt(prediction_error)
    #MAE
    absolute_error <- mean(abs(y.validate-predictions))

    split.robust.lasso.vector <- c(err.curves,bestlambda,prediction_error,root_prediction_error,absolute_error)
    return(split.robust.lasso.vector)



  }

}




