#'Internal lassoenet functions
#'
#'Internal lassoenet functions
#'
#' @param best.lasso.result A vector contains the outputs from the function \code{\link{prediction_Lasso}}.\cr
#' @param best.EN.result A vector contains the outputs from the function \code{\link{prediction_ElasticNet}}.\cr
#' @param data A well-cleaned \code{data.frame}.\cr
#' @param x.indices Locations of the predictors within the \code{data.frame}.\cr
#' @param response Location of the response within the \code{data.frame}.\cr
#' @param alph An internal argument.\cr
#' @param parallel Parallelsation.\cr
#' @return A vector of outputs and some plots related to the best model.\cr
#'
#' @section Details: These are not intended for use by users.
#' This function provides an automatic comparison between the Lasso and the Elastic Net model based on the  Mean squared errors of these models. The return from this function will then enter the \code{\link{penalised_pred}} function.
#'
#'@author Mokyo Zhou
#'
#'
#'@export
#'


prediction.nonsplit.result <- function(best.lasso.result = best.lasso.result, best.EN.result = best.EN.result,
                                       data = data, x.indices = x.indices, response=response, alph = alph, parallel = parallel){

  dataaa <- convertor(data = data, response = response, x.indices = x.indices)
  x <- dataaa[[2]]
  y <- dataaa[[1]]
  #getting the best
  best.lasso.lambda <- best.lasso.result[2]
  best.EN.lambda <- best.EN.result[3]
  best.EN.alpha <- best.EN.result[2]

  #comparision of lasso vs elastic net

  #lasso wins
  if(best.lasso.result[1]==1234567){
    #graphs
    set.seed(1234567); bestlasso.cv <- glmnet::cv.glmnet(x,y, alpha = alph,parallel = parallel)
    bestlasso.glmnet <- glmnet::glmnet(x,y,alpha=alph)
    plot(bestlasso.cv)
    graphics::title("MSE vs log(lambdas)",line=2.2)
    plot(bestlasso.glmnet,xvar="lambda",label = TRUE)
    graphics::title("Regularisation paths of the covariates",line=2.2)
    if(alph == 1){
      abline(v=log(best.lasso.lambda),lty=2)
      a <- glmnet::glmnet(x,y, lambda = best.lasso.lambda)
      mo.formula <- a$call
      mo.coef <- data.frame(as.matrix(coef(a)));rownames(mo.coef)[2:nrow(mo.coef)] <- colnames(x);rownames(mo.coef)[1] <- "intercept"
      colnames(mo.coef) <- "coefficient estimates"
      proportion.explained <- a$dev.ratio
      result.table <- data.frame(best.lasso.result[2], proportion.explained, best.lasso.result[3],
                                 best.lasso.result[4],best.lasso.result[5], best.lasso.result[1])
      rownames(result.table) <- "lasso (non-split, non-repeat)"
      colnames(result.table) <- c("optimum lambda","%null deviance explained","MSE","MSE 95% lower","MSE 95% upper","seed")
      mo.call <- list("Lasso model call" = mo.formula, "coefficients" = mo.coef, "summaries" = result.table)
      return(mo.call)

    }else{
      abline(v=log(best.EN.lambda),lty=2)
      a <- glmnet :: glmnet(x,y, lambda = best.EN.lambda,alpha=best.EN.alpha)
      mo.formula <- a$call
      mo.coef <- data.frame(as.matrix(coef(a)));rownames(mo.coef)[2:nrow(mo.coef)] <- colnames(x);rownames(mo.coef)[1] <- "intercept"
      colnames(mo.coef) <- "coefficient estimates"
      proportion.explained <- a$dev.ratio
      result.table <- data.frame(best.EN.result[2], best.EN.result[3], proportion.explained,
                                 best.EN.result[4],best.EN.result[5],best.EN.result[6], best.EN.result[1])
      rownames(result.table) <- "ElasticNet (non-split, non-repeat)"
      colnames(result.table) <- c("optimum alpha","optimum lambda","%null deviance explained","MSE","MSE 95% lower","MSE 95% upper","seed")
      mo.call <- list("ElasticN model call" = mo.formula, "coefficients" = mo.coef, "summaries" = result.table)
      return(mo.call)
    }

  }else if (best.lasso.result[1] != 1234567){
    bestlasso.glmnet <- glmnet :: glmnet(x,y,alpha=alph)
    if(alph == 1){
      plot(bestlasso.glmnet,xvar="lambda",label = TRUE)
      graphics::title("best model:Lasso",line= 2.2)
      abline(v=log(best.lasso.lambda),lty=2)
      a <- glmnet::glmnet(x,y, lambda = best.lasso.lambda)
      mo.formula <- a$call
      mo.coef <- data.frame(as.matrix(coef(a)));rownames(mo.coef)[2:nrow(mo.coef)] <- colnames(x);rownames(mo.coef)[1] <- "intercept"
      colnames(mo.coef) <- "coefficient estimates"
      proportion.explained <- a$dev.ratio
      result.table <- data.frame(best.lasso.result[2], proportion.explained, best.lasso.result[3],
                                 best.lasso.result[4],best.lasso.result[5], best.lasso.result[1])
      rownames(result.table) <- "lasso (non-split, repeated)"
      colnames(result.table) <- c("optimum lambda","%null deviance explained","MSE","MSE 95% lower","MSE 95% upper","number of error curves")
      mo.call <- list("Lasso model call" = mo.formula, "coefficients" = mo.coef, "summaries" = result.table)
      return(mo.call)

    }else{
      plot(bestlasso.glmnet,xvar="lambda",label = TRUE)
      graphics::title("best model:Elastic Net",line=2.2)
      abline(v=log(best.EN.lambda),lty=2)
      a <- glmnet::glmnet(x,y, lambda = best.EN.lambda,alpha=best.EN.alpha)
      mo.formula <- a$call
      mo.coef <- data.frame(as.matrix(coef(a)));rownames(mo.coef)[2:nrow(mo.coef)] <- colnames(x);rownames(mo.coef)[1] <- "intercept"
      colnames(mo.coef) <- "coefficient estimates"
      proportion.explained <- a$dev.ratio
      result.table <- data.frame(best.EN.result[2], best.EN.result[3], proportion.explained,
                                 best.EN.result[4],best.EN.result[5],best.EN.result[6], best.EN.result[1])
      rownames(result.table) <- "ElasticNet (non-split, repeated)"
      colnames(result.table) <- c("optimum alpha","optimum lambda","%null deviance explained","MSE","MSE 95% lower","MSE 95% upper","number of error curves")
      mo.call <- list("ElasticN model call" = mo.formula, "coefficients" = mo.coef, "summaries" = result.table)
      return(mo.call)

    }

  }


}
