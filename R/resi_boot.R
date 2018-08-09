#' Internal lassoenet function
#'
#' Internal lassoenet function
#'
#'
#'
#' @param x A \code{model.matrix} of the predictors.\cr
#' @param y A vector of response values for the model fitting.\cr
#' @param B.rep The number of residual bootstrappings to do.
#' @param significance The significance level of the confidence intervals e.g. 100(1-\eqn{\alpha})\%.\cr
#' @param best.lambda The optimum lambda selected from either the "OLS" or "ridge" method.
#' @param best.coef A vector of coefficients from either the "OLS" or "ridge" method for constructing inital weights.\cr
#' @param bestgamma The optimum gamma selected from either the "OLS" or "ridge" method.\cr
#' @param ridgelambda The optimum lambda for the ridge regression (If the method "ridge" has been selecteed).\cr
#' @param parallel Parallelisation
#' @param method 1 indicates "OLS", 2 indicates "ridge".
#' @param dataa Your full \code{data.frame}
#' @param xx.indices Locations of the predictors within dataa.
#'
#'
#'
#' @return The 100(1-\eqn{\alpha})\% confidence intervals for the point estimates from using residual bootstrapping.
#'
#' @section Details: These are not intended for use by users. This function takes in outputs from one of \code{\link{OLS.approach}} or \code{\link{OLS.robust}} or \code{\link{ridge.approach}} or \code{\link{ridge.robust}} and
#' computes the corresponding 100(1-\eqn{\alpha})\% confidence intervals for the parameters. The return from this function will enter the funtion \code{\link{Adlasso}} for futher wrapping.
#'
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



resi_boot <- function(x = x, y = y, B.rep = B.rep, significance = significance, best.lambda = best.lambda, best.coef =best.coef,
                      bestgamma=bestgamma, ridgelambda = ridgelambda, parallel = parallel,
                      method = 1, dataa=data, xx.indices = c(1,2,3)){

  #getting fitted values
  fitted <- x %*% best.coef[-1]
  fitted <- fitted + best.coef[1]

  #residuals
  resi <-as.numeric(y)-fitted
  resi <- resi - mean(resi)


  sig <- c(significance/2,1-(significance/2))


  `%dopar%` <- foreach::`%dopar%`

  result.matrix <- foreach::foreach(i = 1:B.rep, .combine = cbind) %dopar%{
    #resampling
    resample <- sample(resi,length(y),T)

    #new response
    new.y <- fitted + resample

    #obtaining new weights
    if (method == 1){
      xl <- dataa[,xx.indices]
      new.data <- data.frame(new.y, xl)
      new.OLS.coef <- summary(lm(new.y~.,data=new.data))$coefficients[,1]

      #new weights
      new.weight <- 1/abs(new.OLS.coef[2:length(new.OLS.coef)])^bestgamma
      new.weight[which(new.weight == Inf)] <- 999999999

      #refit the adpative lasso with new OLS weight
      new.OLS.ad <- glmnet::glmnet(x, new.y, penalty.factor = new.weight, lambda = best.lambda)
      new.adOLS.coef <- coef(new.OLS.ad, s=best.lambda)
      new.adOLS.coef <- as.vector((as.matrix(new.adOLS.coef)))

      data.frame(new.adOLS.coef)
    }else if(method == 2){

      new.ridge <- glmnet::glmnet(x, new.y, alpha = 0, lambda = ridgelambda)
      new.ridge.coef <- coef(new.ridge, s = ridgelambda)

      new.weight <- 1/abs(new.ridge.coef[2:length(new.ridge.coef)])^bestgamma
      new.weight[which(new.weight == Inf)] <- 999999999

      new.ridge.ad <- glmnet::glmnet(x, new.y, penalty.factor = new.weight, lambda = best.lambda)
      new.adridge.coef <- coef(new.ridge.ad, s=best.lambda)
      new.adridge.coef <- as.vector((as.matrix(new.adridge.coef)))

      data.frame(new.adridge.coef)

    }

  }

  result.matrix <- t(result.matrix)
  #CIs construction
  CIs.quantiles <- t(apply(result.matrix, 2, function(x) quantile(x, c(sig[1],sig[2]))))
  CIs.low <- 2 * (best.coef) - CIs.quantiles[,2]; CIs.low <- as.matrix(CIs.low)
  CIs.up <- 2 * (best.coef) - CIs.quantiles[,1];  CIs.up <- as.matrix(CIs.up)
  CI.matrix <- data.frame(CIs.low,CIs.up)
  CI.0.index1 <- which(CI.matrix[,1]==0)
  CI.0.index2 <- which(CI.matrix[,2]==0)
  CI.0.index <- unique(c(CI.0.index1,CI.0.index2))
  zero <- apply(result.matrix,2, function(x) mean(x==0))
  zero <- 1 - zero
  CIs <- data.frame(round(as.vector(best.coef),3),round(as.vector(CIs.low),3), round(as.vector(CIs.up),3),round(as.vector(zero),3))
  if(length(CI.0.index)!=0){
    cat("                                                                                              ")
    CIs[CI.0.index,c(2,3)] <- c(".",".")

  }
  ll <- paste(c(sig[1],sig[2]))
  colnames(CIs) <- c("Estimates",ll[1], ll[2], "proportion of non-zero estimates")
  rownames(CIs)[2:nrow(CIs)] <- colnames(x)
  rownames(CIs)[1] <- "intercept"



  return(CIs)
}



