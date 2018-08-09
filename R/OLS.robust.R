#' Internal lassoenet function
#'
#' Internal lassoenet function
#'
#'
#'
#' @param x A \code{model.matrix} of the predictors.\cr
#' @param y A vector of response values for the model fitting.\cr
#' @param OLS.coef A vector of coefficients from an \code{lm()} object.\cr
#' @param gamma.seq The gamma sequence to try.\cr
#' @param B.rep The number of residual bootstrappings to do.
#' @param significance The significance level of the confidence intervals e.g. 100(1-\eqn{\alpha})\%.\cr
#' @param err.curves The number of error curves to be fitted. Default is 0.\cr
#' @param result.matrix A matrix for storing results.\cr
#' @param CI TRUE for residual bootstrapping. In this version is always TRUE.
#' @param dataa Your full \code{data.frame}
#' @param xx.indices Locations of the predictors within dataa.
#' @param parallel Parallelisation
#'
#' @return A vector of outputs of the best Adpative Lasso model and the 100(1-\eqn{\alpha})\% confidence intervals for the point estimates from using residual bootstrapping.
#'
#' @section Details: These are not intended for use by users. This function is one of the main engines for the Adaptive Lasso computation. This function is used when the weighting method "OLS" is selected.
#' This together with the function \code{\link{OLS.approach}} forms the computation operator for the Adaptive Lasso when using the weighting method "OLS". The return from this function will enter \code{\link{Adlasso}} for futher wrapping.
#'
#' @importFrom stats aggregate coef lm model.matrix predict quantile
#' @importFrom utils data
#' @importFrom graphics abline lines plot points
#'
#' @author Mokyo Zhou
#'
#'
#'
#' @export
#'
#'
#'









OLS.robust <- function(x = x, y = y, OLS.coef = c(1,2,3), gamma.seq = gamma.seq, B.rep = B.rep, significance = significance,
                       err.curves = err.curves, result.matrix = result.matrix, CI = CI, dataa=dataa, xx.indices = xx.indices,
                       parallel = parallel){
  gamma.len <- length(gamma.seq)
  cat("                                                                                               ")


  for(j in 1:gamma.len){
    aa <- gamma.seq[j]

    penalty.vector <- 1/abs(OLS.coef[2:length(OLS.coef)])^gamma.seq[j]

    penalty.vector[which(penalty.vector == Inf)] <- 999999999
    `%dopar%` <- foreach::`%dopar%`




    lambdas <- foreach::foreach(i=1:err.curves, .combine = rbind) %dopar%{
      fit <- glmnet::cv.glmnet(x,y,parallel=parallel,penalty.factor = penalty.vector)
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

  adOLS.robust.vector <- c(err.curves,optimum_pair)

  #reproducing the best coef ########################################################################
  penalty.vector <- 1/abs(OLS.coef[2:length(OLS.coef)])^adOLS.robust.vector[2]                      #
  penalty.vector[which(penalty.vector == Inf)] <- 999999999                                         #
  bestcoeff <- glmnet::glmnet(x, y, penalty.factor = penalty.vector,                   #
                              lambda=adOLS.robust.vector[3])                                                 #
  bestcoef <- coef(bestcoeff)                                                                       #
  proportion.explain <- bestcoeff$dev.ratio                                                         #
  ###################################################################################################

  adOLS.robust.vector <- c(err.curves,optimum_pair,proportion.explain)
  adOLS.robust.vector <- base::unname(adOLS.robust.vector)
  if(CI == TRUE) {


    CIs <- resi_boot(x = x, y = y, best.lambda = adOLS.robust.vector[3], best.coef =bestcoef,
                     bestgamma=adOLS.robust.vector[2], ridgelambda = NULL, parallel = parallel, method = 1, dataa=dataa, xx.indices = xx.indices, B.rep= B.rep, significance = significance)


  }



  return(list(adOLS.robust.vector,CIs))
}



