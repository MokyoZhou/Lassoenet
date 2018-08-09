#' Internal lassoenet function
#'
#' Internal lassoenet function
#'
#' @param err.curves The number of error curves to be fitted. Default is 0.\cr
#' @param x A model.matrix for the predictors.\cr
#' @param y A vector of response values.\cr
#' @param type.lambda Either "lambda.min" or "lambda.1se", default is "lambda.min.\cr
#' @param al 1 for the Lasso and 0 for ridge regression.\cr
#' @param parallel Parallelisation
#' @return A vector of results for the best Lasso model under the condition where the full dataset has been used for modelling. The return from this function will enter \code{\link{prediction_Lasso}}.
#'
#' @section Details: These are not intended for use by users. This function is one of the main engines for the Lasso computation. This function is used when the user does not want to split the dataset into
#'  a training and a test set. This together with the function \code{\link{split_lasso}} form the overall computation operator for the Lasso models. The return from this function will enter \code{\link{prediction_Lasso}} for futher wrapping.
#'
#'
#'
#' @author Mokyo Zhou
#'
#'
#'
#' @export
#'
no_split_lasso <- function(err.curves = 0, x = x, y=y,type.lambda=type.lambda, al = 1, parallel = parallel){

  if(err.curves==0){

    z <- ifelse(type.lambda=="lambda.min",z <- 9, z <- 10)

    set.seed(1234567);best.lambda<-glmnet::cv.glmnet(x, y, parallel = parallel )
    prediction_error <- best.lambda$cvm[which(best.lambda$lambda==best.lambda[z])]
    prediction_up <- prediction_error + 1.96 * best.lambda$cvsd[which(best.lambda$lambda==best.lambda[z])]
    prediction_low <- prediction_error - 1.96 * best.lambda$cvsd[which(best.lambda$lambda==best.lambda[z])]

    ifelse(z==9,bestLambda <- best.lambda$lambda.min, bestLambda <- best.lambda$lambda.1se)
    nosplit.lasso.vector <- c(1234567, bestLambda,prediction_error,prediction_low,prediction_up)
    nosplit.lasso.vector <- base::unname(nosplit.lasso.vector)
    return(nosplit.lasso.vector)

  }else{
    cat("                                                                                               ")

    `%dopar%` <- foreach::`%dopar%`
    lambdas.full <- foreach::foreach(i=1:err.curves, .combine = rbind) %dopar%{
      fit <- glmnet::cv.glmnet(x,y,alpha=al,parallel=parallel)
      data.frame(fit$lambda, fit$cvm)
    }



    # take mean cvm for each lambda
    lambdas <- lambdas.full
    lambdas <- aggregate(lambdas[, 2], list(lambdas$fit.lambda), mean)
    if(al!=0){
      plot(log(lambdas.full[,1]),lambdas.full[,2],pch=20,col="black",main="Averaging the error curves (Lasso)",xlab="log lambdas",
           ylab="MSE")
      lines(log(lambdas[,1]),lambdas[,2],col="blue")
    }

    # select the best one
    bestindex <- which(lambdas[2]==min(lambdas[2]))
    bestlambda <- lambdas[bestindex,1]
    prediction_error <- lambdas[bestindex,2]

    best_lambda_indices <- which(lambdas.full[1]==bestlambda)
    best_cvms <- lambdas.full[best_lambda_indices,2]
    prediction_CI <- quantile(best_cvms,c(0.025,0.975))

    nosplit.lasso.robust.vector <- c(err.curves,bestlambda,prediction_error,prediction_CI)
    nosplit.lasso.robust.vector <- base::unname(nosplit.lasso.robust.vector)
    if(al!=0){
      abline(v=log(nosplit.lasso.robust.vector[2]),lty=2)
      a <- (rep(log(nosplit.lasso.robust.vector[2]),3)); b<- nosplit.lasso.robust.vector[c(3,4,5)];points(a,b,col="red",pch=20)

    }


    return(nosplit.lasso.robust.vector)

  }

}
