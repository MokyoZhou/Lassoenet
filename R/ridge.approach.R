#' Internal lassoenet function
#'
#' Internal lassoenet function
#'
#'
#'
#' @param x A \code{model.matrix} of the predictors.\cr
#' @param y A vector of response values for the model fitting.\cr
#' @param err.curves The number of error curves to be fitted. Default is 0.\cr
#' @param type.lambda Either "lambda.min" or "lambda.1se".
#' @param gamma.seq The gamma sequence to try.\cr
#' @param CI TRUE for residual bootstrapping. In this version is always TRUE.
#' @param dataa Your full \code{data.frame}
#' @param parallel Parallelisation
#' @param xx.indices Locations of the predictors within dataa.
#' @param B.rep The number of residual bootstrappings to do.
#' @param significance The significance level of the confidence intervals e.g. 100(1-\eqn{\alpha})\%.\cr
#' @return A vector of outputs of the best Adpative Lasso model and the 100(1-\eqn{\alpha})\% confidence intervals for the point estimates from using residual bootstrapping.
#'
#' @section Details: These are not intended for use by users. This function is one of the main engines for the Adaptive Lasso computation. This function is used when the weighting method "ridge" is selected.
#' This together with the function \code{\link{ridge.robust}} forms the computation operator for the Adaptive Lasso when using the weighting method "ridge". The return from this function will enter \code{\link{Adlasso}} for futher wrapping.
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






ridge.approach <- function(x = x, y = y, err.curves = err.curves, type.lambda = type.lambda, gamma.seq = c(0.5,1,2), CI = CI,
                          dataa = dataa, parallel = parallel, xx.indices = xx.indices, B.rep = B.rep, significance = significance){

  i <- NULL

  if(err.curves==0){

    #ridge coefficients
    set.seed(1234567)
    ridge.reg <- glmnet::cv.glmnet(x,y, alpha=0,parallel = parallel)

    z <- ifelse(type.lambda=="lambda.min",z <- 9, z <- 10)

    ridge.coefficients <- coef(ridge.reg, s=ridge.reg$lambda[which(ridge.reg$lambda==ridge.reg[z])])
    ridgelambda <- ridge.reg$lambda[which(ridge.reg$lambda==ridge.reg[z])]
    `%dopar%` <- foreach::`%dopar%`




    result.matrix <-foreach::foreach(i = gamma.seq, .combine = rbind) %dopar% {

      penalty.vector <- 1/abs(ridge.coefficients[2:length(ridge.coefficients)])^i


      penalty.vector[which(penalty.vector == Inf)] <- 999999999

      set.seed(1234567);cv <-glmnet::cv.glmnet(x,y, parallel = parallel, penalty.factor = penalty.vector)
      prediction_error <- cv$cvm[which(cv$lambda==cv[z])]
      prediction_up <- prediction_error + 1.96 * cv$cvsd[which(cv$lambda==cv[z])]
      prediction_low <- prediction_error - 1.96 * cv$cvsd[which(cv$lambda==cv[z])]
      data.frame(i, cv[z], prediction_error,prediction_low,prediction_up)
    }
    result.matrix <- as.matrix(result.matrix)


    optimum_pair <-result.matrix[which(result.matrix[,3]== min(result.matrix[,3])),]

    ad.nonerobust.vector <- c(1234567, optimum_pair)

    #reproducing the best coef ########################################################################
    penalty.vector <- 1/abs(ridge.coefficients[2:length(ridge.coefficients)])^ad.nonerobust.vector[2] #
    penalty.vector[which(penalty.vector == Inf)] <- 999999999                                         #
    bestcoeff <- glmnet::glmnet(x, y, penalty.factor = penalty.vector,                   #
                       lambda=ad.nonerobust.vector[3])                                                #
    bestcoef <- coef(bestcoeff)                                                                       #
    proportion.explain <- bestcoeff$dev.ratio                                                         #
    ###################################################################################################

    ad.nonerobust.vector <- c(1234567, optimum_pair,ridgelambda,proportion.explain)
    ad.nonerobust.vector <- base::unname(ad.nonerobust.vector)

    if(CI == TRUE){

      CIs <- resi_boot(x = x, y = y, best.lambda = ad.nonerobust.vector[3], best.coef =bestcoef,
                       bestgamma=ad.nonerobust.vector[2], ridgelambda = ridgelambda , parallel = parallel,
                       method = 2, dataa = dataa, xx.indices = xx.indices, B.rep = B.rep, significance = significance)


    }
    return(list(ad.nonerobust.vector,CIs))
  }else{
    result.matrix <- matrix(0,length(gamma.seq),5)
    robust.ridge <- ridge.robust(x = x, y = y, type.lambda=type.lambda, al=0, gamma.seq = gamma.seq, err.curves = err.curves,
                                 result.matrix = result.matrix, CI=CI,
                                 dataa = dataa, xx.indices = xx.indices, parallel = parallel, B.rep=B.rep, significance=significance)
    return(robust.ridge)

  }

}
