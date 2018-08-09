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
#' @param splits A two element vector with the first element being the training proportion and the second element being the testing proportion.\cr
#' @param parallel Parallelsation.\cr
#' @return A vector of outputs and some plots related to the best model.
#'
#' @section Details: These are not intended for use by users.
#' This function provides an automatic comparison between the Lasso and the Elastic Net model based on the user's preferred prediction metric. The return from this function will then enter the \code{\link{penalised_pred}} function.\cr
#'
#'@author Mokyo Zhou
#'
#'
#'@export
#'










prediction.split.result <- function(best.lasso.result = best.lasso.result, best.EN.result = best.EN.result,
                                    data = data, x.indices = x.indices, response=response, alph = alph,
                                    splits = splits,parallel = parallel){


  split.data <- datasplit(data = data, splits = splits)
  dataaa <- convertor(data = split.data[[1]], response = response, x.indices = x.indices)
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
    set.seed(1234567); bestlasso.cv <- glmnet::cv.glmnet(x,y, alpha = alph,parallel = parallel)#train
    bestlasso.glmnet <- glmnet::glmnet(x,y,alpha=alph)#train

    plot(bestlasso.cv)
    graphics::title("MSE vs log.lambda (on train)", line =2.2)
    plot(bestlasso.glmnet,xvar="lambda",label = TRUE)
    graphics::title("Regularisation paths of the covariates (on train)",line = 2.2)

    if(alph == 1){
      abline(v=log(best.lasso.lambda),lty=2)
      a <- glmnet::glmnet(x,y, lambda = best.lasso.lambda)
      mo.formula <- a$call
      mo.coef <- data.frame(as.matrix(coef(a)));rownames(mo.coef)[2:nrow(mo.coef)] <- colnames(x);rownames(mo.coef)[1] <- "intercept"
      colnames(mo.coef) <- "coefficient estimates (on train)"
      proportion.explained <- a$dev.ratio
      result.table <- data.frame(best.lasso.result[2], proportion.explained, best.lasso.result[3],
                                 best.lasso.result[4],best.lasso.result[5], best.lasso.result[1])
      rownames(result.table) <- "lasso (split, non-repeat) "
      colnames(result.table) <- c("optimum lambda","%null deviance explained (On Train)","out of sample MSE",
                                  "out of sample RMSE","out of sample MAE","seed for both data splitting and the glmnet fit")
      mo.call <- list("Lasso model call" = mo.formula, "coefficients" = mo.coef, "summaries" = result.table)
      return(mo.call)
    }else{
      abline(v=log(best.EN.lambda),lty=2)
      a <- glmnet::glmnet(x,y, lambda = best.EN.lambda,alpha=best.EN.alpha)
      mo.formula <- a$call
      mo.coef <- data.frame(as.matrix(coef(a)));rownames(mo.coef)[2:nrow(mo.coef)] <- colnames(x);rownames(mo.coef)[1] <- "intercept"
      colnames(mo.coef) <- "coefficient estimates (on train)"
      proportion.explained <- a$dev.ratio
      result.table <- data.frame(best.EN.result[2], best.EN.result[3], proportion.explained,
                                 best.EN.result[4],best.EN.result[5],best.EN.result[6], best.EN.result[1])
      rownames(result.table) <- "ElasticNet (split, non-repeat)"
      colnames(result.table) <- c("optimum alpha","optimum lambda","%null deviance explained (On Train)",
                                  "out of sample MSE","out of sample RMSE","out of sample MAE","seed for both data splitting and the glmnet fit")
      mo.call <- list("ElasticN model call" = mo.formula, "coefficients" = mo.coef, "summaries" = result.table)
      return(mo.call)
    }


  }else if (best.lasso.result[1] != 1234567){
    bestlasso.glmnet <- glmnet::glmnet(x,y,alpha=alph)
    if(alph == 1){
      plot(bestlasso.glmnet,xvar="lambda",label = TRUE)
      graphics::title("best model:Lasso", line=2.2)
      abline(v=log(best.lasso.lambda),lty=2)
      a <- glmnet::glmnet(x,y, lambda = best.lasso.lambda)
      mo.formula <- a$call
      mo.coef <- data.frame(as.matrix(coef(a)));rownames(mo.coef)[2:nrow(mo.coef)] <- colnames(x);rownames(mo.coef)[1] <- "intercept"
      colnames(mo.coef) <- "coefficient estimates (on train)"
      proportion.explained <- a$dev.ratio
      result.table <- data.frame(best.lasso.result[2], proportion.explained, best.lasso.result[3],
                                 best.lasso.result[4],best.lasso.result[5], best.lasso.result[1])
      rownames(result.table) <- "lasso (split, repeated)"
      colnames(result.table) <- c("optimum lambda","%null deviance explained (On Train)",
                                  "out of sample MSE","out of sample RMSE","out of sample MAE","number of error curves")
      mo.call <- list("Lasso model call" = mo.formula, "coefficients" = mo.coef, "summaries" = result.table)
      return(mo.call)
    }else{
      plot(bestlasso.glmnet,xvar="lambda",label = TRUE)
      graphics::title("best model:Elastic Net", line= 2.2)
      abline(v=log(best.EN.lambda),lty=2)
      a <- glmnet::glmnet(x,y, lambda = best.EN.lambda,alpha=best.EN.alpha)
      mo.formula <- a$call
      mo.coef <- data.frame(as.matrix(coef(a)));rownames(mo.coef)[2:nrow(mo.coef)] <- colnames(x);rownames(mo.coef)[1] <- "intercept"
      colnames(mo.coef) <- "coefficient estimates (on train)"
      proportion.explained <- a$dev.ratio
      result.table <- data.frame(best.EN.result[2], best.EN.result[3], proportion.explained,
                                 best.EN.result[4],best.EN.result[5],best.EN.result[6], best.EN.result[1])
      rownames(result.table) <- "ElasticNet (split, repeated)"
      colnames(result.table) <- c("optimum alpha","optimum lambda","%null deviance explained (On Train)",
                                  "out of sample MSE","out of sample RMSE","out of sample MAE","number of error curves")
      mo.call <- list("ElasticN model call" = mo.formula, "coefficients" = mo.coef, "summaries" = result.table)
      return(mo.call)
    }

  }




}




