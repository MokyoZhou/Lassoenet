#' Internal lassoenet function
#'
#' Internal lassoenet function
#'
#' @param x A \code{model.matrix} of the preditors for the training data.\cr
#' @param y A vector of response values on the training set.\cr
#' @param test A \code{model.matrix} of the predictors for the testing data.\cr
#' @param yv A vector of response values on the testing set.\cr
#' @param alpha.seq The alpha sequence for the Elastic Net.\cr
#' @param err.curves The number of error curves to be fitted. Default is 0.\cr
#' @param result.matrix A matrix for storing results.\cr
#' @param parallel Parallelisation
#' @return A vector of results for the best Elastic Net model along with some plots of the error curves, under the condition where the full dataset has been divided into a traning and a test set. The return from this function will enter \code{\link{prediction_ElasticNet}}.
#'
#' @section Details: These are not intended for use by users. This function is one of the main engines for the Elastic Net computation. This function is used when the user would like to split the dataset into
#'  a training and a test set. This together with the function \code{\link{split_en}} forms the computation operator for the Elastic Net when using the a traning and a testing set. The return from this function will enter \code{\link{prediction_ElasticNet}} for futher wrapping.
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

en.robust.split <- function(x = x, y = y, test = test, yv =yv, alpha.seq = alpha.seq, err.curves = err.curves,
                            result.matrix = result.matrix, parallel = parallel){
  alpha.len <- length(alpha.seq)
  full.matrix <- NULL
  reduced.matrix <- NULL
  `%dopar%` <- foreach::`%dopar%`

  #breaking the inherent seed within `datasplit`
  set.seed(Sys.time())
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
    bestindex <- which(lambdas[2]==min(lambdas[2]))
    bestlambda <- lambdas[bestindex,1]
    prediction_error <- lambdas[bestindex,2]


    result.matrix[j,] <- unlist(c(aa, bestlambda, prediction_error))



  }
  optimum_index <- which(result.matrix[,3]== min(result.matrix[,3]))
  optimum_pair <-result.matrix[optimum_index,]

  full.vec <- which(full.matrix[,1]==0);full.vec <- c(1,full.vec)
  reduced.vec <- which(reduced.matrix[,1]==0);reduced.vec <- c(1, reduced.vec)

  best.full <- full.matrix[full.vec[optimum_index]:full.vec[optimum_index+1],]
  zeroind <- which(best.full[,1]==0); best.full <- best.full[-zeroind,]
  best.reduce <- reduced.matrix[reduced.vec[optimum_index]:reduced.vec[optimum_index+1],]
  zeroind <- which(best.reduce[,1]==0); best.reduce <- best.reduce[-zeroind,]


  plot(log(best.full[,1]),best.full[,2],pch=20,col="black",main="Averaging across the error curves (EN, Training Set)"
       ,xlab="log lambdas",ylab="MSE")
  lines(log(best.reduce[,1]),best.reduce[,2],col="green")
  abline(v=log(optimum_pair[2]), lty=2)


  best.en <- glmnet::glmnet(x,y,alpha=optimum_pair[1])

  predictions <- predict(best.en, newx = test, s = optimum_pair[2])



  #MSE
  prediction_error <- sum((yv - predictions)^2)/nrow(test)
  #Root MSE
  root_prediction_error <- sqrt(prediction_error)
  #MAE
  absolute_error <- mean(abs(yv-predictions))


  split.en.robust.vector <- c(err.curves,optimum_pair[1],optimum_pair[2],prediction_error,root_prediction_error,absolute_error)
  split.en.robust.vector <- base::unname(split.en.robust.vector)
  return(split.en.robust.vector)


}




