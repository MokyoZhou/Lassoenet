#' Internal lassoenet function
#'
#' Internal lassoenet function
#'
#'
#'
#' @param x A \code{model.matrix} of the predictors.\cr
#' @param y A vector of response values for the model fitting.\cr
#' @param type.lambda Either "lambda.min" or "lambda.1se".\cr
#' @param al 0 for ridge.\cr
#' @param gamma.seq The gamma sequence to try.\cr
#' @param err.curves The number of error curves to be fitted. Default is 0.\cr
#' @param result.matrix A matrix for storing the results.\cr
#' @param CI TRUE for residual bootstrapping. In this version is always TRUE.
#' @param B.rep The number of residual bootstrappings to do.
#' @param significance The significance level of the confidence intervals e.g. 100(1-\eqn{\alpha})\%.\cr
#' @param dataa Your full \code{data.frame}
#' @param xx.indices Locations of the predictors within dataa.
#' @param parallel Parallelisation
#' @return A vector of outputs of the best Adpative Lasso model and the 100(1-\eqn{\alpha})\% confidence intervals for the point estimates from using residual bootstrapping.
#'
#' @section Details: These are not intended for use by users. This function is one of the main engines for the Adaptive Lasso computation. This function is used when the weighting method "ridge" is selected.
#' This together with the function \code{\link{ridge.approach}} forms the computation operator for the Adaptive Lasso when using the weighting method "ridge". The return from this function will enter \code{\link{Adlasso}} for futher wrapping.
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




ridge.robust <- function(x = x, y = y, type.lambda=type.lambda, al=0, gamma.seq = gamma.seq, err.curves = err.curves,
                         result.matrix = result.matrix, CI=CI, B.rep = B.rep, significance = significance,
                         dataa = dataa, xx.indices = xx.indices, parallel = parallel){

  gamma.len <- length(gamma.seq)




  #ridge loop
  ridge.lambda <- no_split_lasso(x = x, y = y, err.curves = err.curves, type.lambda=type.lambda, al = al,parallel = parallel)
  ridge.lambda <- ridge.lambda[2]


  #ridge coef
  ridge.coef <- coef(glmnet::glmnet(x,y, alpha=0, lambda = ridge.lambda))
  `%dopar%` <- foreach::`%dopar%`

  for(j in 1:gamma.len){
    aa <- gamma.seq[j]

    penalty.vector <- 1/abs(ridge.coef[2:length(ridge.coef)])^gamma.seq[j]

    penalty.vector[which(penalty.vector == Inf)] <- 999999999



    lambdas <- foreach::foreach(i=1:err.curves, .combine = rbind) %dopar%{
      fit <- glmnet::cv.glmnet(x,y,penalty.factor = penalty.vector,parallel=parallel)
      data.frame(fit$lambda, fit$cvm)
    }

    # take mean cvm for each lambda
    lambdas.full <- lambdas
    lambdas <- aggregate(lambdas[, 2], list(lambdas$fit.lambda), mean)

    # select the best one
    bestindex <- which(lambdas[2]==min(lambdas[2]))
    bestlambda <- lambdas[bestindex,1]
    prediction_error <- lambdas[bestindex,2]

    best_lambda_indices <- which(lambdas.full[,1]==bestlambda)
    best_cvms <- lambdas.full[best_lambda_indices,2]
    prediction_CI <- quantile(best_cvms,c(0.025,0.975))

    result.matrix[j,] <- unlist(c(aa, bestlambda, prediction_error, prediction_CI[1], prediction_CI[2]))


  }
  optimum_pair <-result.matrix[which(result.matrix[,3]== min(result.matrix[,3])),]

  adridge.robust.vector <- c(err.curves,optimum_pair,ridge.lambda)

  #reproducing the best coef ########################################################################
  penalty.vector <- 1/abs(ridge.coef[2:length(ridge.coef)])^adridge.robust.vector[2]                #
  penalty.vector[which(penalty.vector == Inf)] <- 999999999                                         #
  bestcoeff <- glmnet::glmnet(x, y, penalty.factor = penalty.vector,                   #
                              lambda=adridge.robust.vector[3])                                               #
  bestcoef <- coef(bestcoeff)                                                                       #
  proportion.explain <- bestcoeff$dev.ratio                                                         #
  ###################################################################################################

  adridge.robust.vector <- c(err.curves,optimum_pair,ridge.lambda,proportion.explain)
  adridge.robust.vector <- base::unname(adridge.robust.vector)

  if(CI == TRUE){


    CIs <- resi_boot(x = x, y = y, best.lambda = adridge.robust.vector[3], best.coef =bestcoef,
                     bestgamma=adridge.robust.vector[2], ridgelambda = ridge.lambda , parallel = parallel, method = 2,
                     dataa = dataa, xx.indices = xx.indices, B.rep= B.rep, significance = significance)


  }



  return(list(adridge.robust.vector,CIs))
}



