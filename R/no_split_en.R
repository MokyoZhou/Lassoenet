#' Internal lassoenet function
#'
#' Internal lassoenet function
#'
#' @param err.curves The number of error curves to be fitted. Default is 0.\cr
#' @param x A model.matrix for the predictors.\cr
#' @param y A vector of response values.\cr
#' @param step.size Step size of the alpha grid.\cr
#' @param type.lambda Either "lambda.min" or "lambda.1se", default is "lambda.min.\cr
#' @param parallel Parallelisation
#' @return A vector of results for the best Elastic Net model under the condition where the full dataset has been used for modelling. The return from this function will enter \code{\link{prediction_ElasticNet}}.
#'
#' @section Details: These are not intended for use by users. This function is one of the main engines for the Elastic Net computation. This function is used when the user does not want to split the dataset into
#'  a training and a test set. This together with the function \code{\link{en.robust}} form the computation operator for the Elastic Net when using the full dataset. The return from this function will enter \code{\link{prediction_ElasticNet}} for futher wrapping.
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







no_split_en <- function(err.curves = 0, x = x, y=y, step.size = 0, type.lambda=type.lambda,parallel = parallel){
  alpha.seq <-seq(0,1,step.size)
  ifelse(alpha.seq[length(alpha.seq)]==1, alpha.seq <- alpha.seq[c(-1,-length(alpha.seq))],
         alpha.seq <- alpha.seq[c(-1)])

  i <- NULL

  if(err.curves==0){
    type.lambda <- type.lambda


    z <- ifelse(type.lambda=="lambda.min",z <- 9, z <- 10)
    `%dopar%` <- foreach::`%dopar%`


    result.matrix <-foreach::foreach(i = alpha.seq, .combine = rbind) %dopar% {

      set.seed(1234567);cv <-glmnet::cv.glmnet(x,y, parallel = parallel, alpha = i)
      prediction_error <- cv$cvm[which(cv$lambda==cv[z])]
      prediction_up <- prediction_error + 1.96 * cv$cvsd[which(cv$lambda==cv[z])]
      prediction_low <- prediction_error - 1.96 * cv$cvsd[which(cv$lambda==cv[z])]
      data.frame(alpha=i, lambda=cv[z], prediction_error,prediction_low,prediction_up)
    }
    result.matrix <- as.matrix(result.matrix)

    optimum_pair <-result.matrix[which(result.matrix[,3]== min(result.matrix[,3])),]

    nosplit.en.vector <- c(1234567, optimum_pair)
    nosplit.en.vector <- base::unname(nosplit.en.vector)
    return(nosplit.en.vector)


  }else if(err.curves!=0){
    result.matrix <-matrix(0,length(alpha.seq),5)
    en.robust.result <- en.robust(x = x, y = y, alpha.seq = alpha.seq, err.curves = err.curves, result.matrix = result.matrix, parallel = parallel)
    return(en.robust.result)

  }

}




