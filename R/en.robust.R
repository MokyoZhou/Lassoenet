#' Internal lassoenet function
#'
#' Internal lassoenet function
#'
#' @param x A model.matrix for the predictors.\cr
#' @param y A vector of response values.\cr
#' @param alpha.seq The alpha sequence for the Elastic Net
#' @param err.curves The number of error curves to be fitted. Default is 0.\cr
#' @param result.matrix A matrix for storing results.\cr
#' @param parallel Parallelisation
#' @return A vector of results for the best Elastic Net model along with some plots of the error curves, under the condition where the full dataset has been used for modelling. The return from this function will enter \code{\link{prediction_ElasticNet}}.
#'
#' @section Details: These are not intended for use by users. This function is one of the main engines for the Elastic Net computation. This function is used when the user does not want to split the dataset into
#'  a training and a test set. This together with the function \code{\link{no_split_en}} form the computation operator for the Elastic Net when using the full dataset. The return from this function will enter \code{\link{prediction_ElasticNet}} for futher wrapping.
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












en.robust <- function(x = x, y = y, alpha.seq = alpha.seq, err.curves = err.curves, result.matrix = result.matrix, parallel = parallel){

  alpha.len <- length(alpha.seq)
  full.matrix <- NULL
  reduced.matrix <- NULL
  `%dopar%` <- foreach::`%dopar%`

  for(j in 1:alpha.len){
    aa <- alpha.seq[j]


    lambdas <- foreach::foreach(i=1:err.curves, .combine = rbind) %dopar%{
      fit <- glmnet::cv.glmnet(x,y,alpha=aa,parallel=parallel)
      data.frame(fit$lambda, fit$cvm)
    }


    # take mean cvm for each lambda
    lambdas.full <- lambdas
    f.matrix <- lambdas.full
    full.matrix <- rbind(full.matrix,f.matrix,c(0,0))

    lambdas <- aggregate(lambdas[, 2], list(lambdas$fit.lambda), mean)
    r.matrix <- lambdas
    reduced.matrix <- rbind(reduced.matrix,r.matrix, c(0,0))


    # select the best one
    bestindex = which(lambdas[2]==min(lambdas[2]))
    bestlambda = lambdas[bestindex,1]
    prediction_error <- lambdas[bestindex,2]

    best_lambda_indices <- which(lambdas.full[,1]==bestlambda)
    best_cvms <- lambdas.full[best_lambda_indices,2]
    prediction_CI <- quantile(best_cvms,c(0.025,0.975))

    result.matrix[j,] <- unlist(c(aa, bestlambda, prediction_error, prediction_CI[1], prediction_CI[2]))


  }
  optimum_index <- which(result.matrix[,3]== min(result.matrix[,3]))
  optimum_pair <-result.matrix[optimum_index,]

  full.vec <- which(full.matrix[,1]==0); full.vec <- c(1, full.vec)
  reduced.vec <- which(reduced.matrix[,1]==0); reduced.vec <- c(1, reduced.vec)

  best.full <- full.matrix[full.vec[optimum_index]:full.vec[optimum_index+1],]
  zeroind <- which(best.full[,1]==0); best.full <- best.full[-zeroind,]
  best.reduce <- reduced.matrix[reduced.vec[optimum_index]:reduced.vec[optimum_index+1],]
  zeroind <- which(best.reduce[,1]==0); best.reduce <- best.reduce[-zeroind,]


  plot(log(best.full[,1]),best.full[,2],pch=20,col="black",main="Averaging the error curves (EN)",
       xlab="log lambdas",ylab="MSE")
  lines(log(best.reduce[,1]),best.reduce[,2],col="green")
  abline(v=log(optimum_pair[2]), lty=2)
  a <- (rep(log(optimum_pair[2]),3)); b<- optimum_pair[c(3,4,5)];points(a,b,col="red",pch=20)

  nosplit.en.robust.vector <- c(err.curves,optimum_pair)
  nosplit.en.robust.vector <- base::unname(nosplit.en.robust.vector)
  return(nosplit.en.robust.vector)

}


