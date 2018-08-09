#'Simulation studies for Lasso, Elastic Net and ridge
#'
#'This function calls the data generating function \code{simulation.generation.data} and can be used
#'to conduct simulation studies for comparison between the Lasso, the Elastic Net and the ridge regression.
#'
#' @param n.resample The number of simulation datasets to generate\cr
#' @param n The number of rows in the each \code{data.frame}
#' @param coeff A vector a true coefficients fixed acrossed the simulated datasets\cr
#' @param matrix.option 1: Using an Exchangeable correlation matrix to simulate the predictors\cr
#'                      2: Using an Autoregressive correlation matrix to simulate the predictors\cr
#' @param collinear The correlation levels within the \code{matrix.option}
#' @param sig The model inherent error, the \eqn{\sigma^2}
#' @param split.prop An element specifying the training proportion. Note the testing proportion will be 1 - the training proportion.
#' @param step.alpha The step size of the alpha grid for the Elastic Net.
#' @param option 1: split the dataset according to \code{c(split.prop, 1 - split.prop)}\cr
#'               2: Use the whole dataset. Note When \code{option = 2}, the \code{split.prop} will be ignored\cr
#'               Note: When using \code{simulation.collinear} please make sure "option = 1" under the current package.
#' @param parallel Parallelisation
#' @return A list of elements:
#' \item{final.table}{A summary table contains:
#'                                              1.The averaged external simulation prediction error for Lasso, Elastic Net and the Ridge\cr
#'                                              2.The standard error of the prediction error by using bootstrapping\cr
#'                                              3.The proportion of times the each method correctly chooses the true model\cr
#'                                              4.The averaged \eqn{\alpha}: 1 for Lasso, averaged simluated \eqn{\alpha} for the Elastic Net and 0 for the Ridge\cr}
#' \item{lasso.result}{A matrix full with simulation results for the Lasso method. Users can freely use this result matrix to obtain further insight\cr}
#' \item{EN.result}{A matrix full with simualtion results for the Elastic Net method. Users can freely use this result matrix to obtain further insight\cr}
#' \item{ridge.result}{A matrix full with simulation results for the ridge regression. Users can freely use this result matrix to botain further isnight\cr}
#'
#' @section Details: The function is one of the core function for the simulation studies. The function supports comparison between the Lasso, the Elastic Net and the ridge regression.
#' This function calls the function \code{simulation.generation.data} and thus, users can study different datasets of their liking. The function provide a summary table for the simulation results.
#' The matrices containing the simulation iterations for the three methods are also provided. Therefore, users are free to conduct further investigation.
#'
#' @examples
#' #number of simulated dataset 20, 200 rows, vector of true coefficient is c(10,8,0,0,12,0,0,0,0,0),
#' #using autoregressive correlation matrix, correlation level is 0.2 in the autoregressive matrix,
#' #model error is 3, training proportion is 0.6, step size of the alpha grid is 0.2, splitting
#' #the dataset into a training and a testing set. No parallelisation.
#' simulation1 <- simulation.collinear(n.resample = 20, n = 200, coeff = c(10,8,0,0,12,0,0,0,0,0),
#' matrix.option = 2, collinear = 0.2, sig = 3,split.prop = 0.6, parallel = FALSE, step.alpha = 0.2,
#' option= 1)
#' @author Mokyo Zhou
#'
#' @export
#'
#'





simulation.collinear <- function(n.resample = n.resample, n = n, coeff = coeff, matrix.option = 2, collinear = collinear, sig = sig,
                                 split.prop = split.prop, step.alpha = step.alpha,option=1,parallel = FALSE){



  Ans <- suppressWarnings(ifelse(is.na(as.numeric(n.resample)),stop("please make sure the number of iterations is a number"),NUL<-0))
  Ans <- ifelse(length(n.resample) !=1,stop("Please make sure the number of iteratiosn is a single number"),NUL<-0)
  Ans <- ifelse(n.resample < 0,stop("please make sure the number of iterations is a positive"),NUL<-0)
  Ans <- ifelse(n.resample != floor(n.resample),stop("please make sure the number of iteration is a whole number"),NUL<-0)
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(n)),stop("please make sure the number of observation is a number"),NUL<-0))
  Ans <- ifelse(length(n) !=1,stop("Please make sure the number of observations is a single number"),NUL<-0)
  Ans <- ifelse(n < 0,stop("please make sure the number of observation is a positive"),NUL<-0)
  Ans <- ifelse(n != floor(n),stop("please make sure the number of observation is a whole number"),NUL<-0)
  Ans <- ifelse(dim(n)[1] < 5 * length(coeff),stop("please make sure you dataset has observations at least 5 times the length of your coeff"),NUL<-0)
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(coeff)),stop("please make sure there are all numbers within coeff"),NUL<-0))
  Ans <- ifelse(length(coeff) ==1,stop("Please make sure you have more than 1 element within the coeff"),NUL<-0)
  Ans <- ifelse(sum(is.na(as.numeric(coeff)))!=0,stop("Please make sure all the elements within coeff are numbers"),NUL<-0)
  Ans <- ifelse(is.na(as.numeric(collinear)),stop("please make sure the collinear is a number"),NUL<-0)
  Ans <- ifelse(length(collinear) != 1,stop("Please make sure you have only entered a single number for collienar"),NUL<-0)
  Ans <- ifelse(collinear < -1 || collinear > 1,stop("Please make sure collinear is within [-1,1]"),NUL<-0)
  Ans <- ifelse(is.na(as.numeric(sig)),stop("please make sure the sig is a number"),NUL<-0)
  Ans <- ifelse(length(sig) != 1,stop("Please make sure you have only entered a single number for sig"),NUL<-0)
  Ans <- ifelse(sig < 0 ,stop("Please make sure the sig is a positive number"),NUL<-0)
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(split.prop)),stop("please only enter a number for the training proportion"),NUL<-0))
  Ans <- ifelse(length(split.prop)!=1,stop("please only enter a number for the training proportion"),NUL<-0)
  Ans <- ifelse(split.prop < 0 || split.prop >1, stop("please make sure the training proportion is within [0,1]"),NUL<-0)
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(step.alpha)),stop("please only enter a number for the step size of alpha"),NUL<-0))
  Ans <- ifelse(length(step.alpha) != 1,stop("please make sure the step size for alpha is a single number"),NUL<-0)
  Ans <- ifelse(step.alpha < 0 || step.alpha >1,stop("please make sure the alpha step size is within [0,1]"),NUL<-0)
  Ans <- ifelse(length(option)!=1,stop("please make sure you have only entered 1 value for option"),NUL<-0)
  Ans <- ifelse(option != 1 & option != 2,stop("Please make sure you have the correct option input, 1:split, 2: nonsplit"),NUL<-0)
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(matrix.option)),stop("please make sure the matrix.option is a number"),NUL<-0))
  Ans <- ifelse(matrix.option !=1 & matrix.option != 2,stop("Please make sure you have the correct matrix.option input, 1:exchangeable, 2:autoregressive"),NUL<-0)
  Ans <- ifelse(length(matrix.option)!=1,stop("Please make sure you have only entered 1 number for the matrix.option input, 1:exchangeable, 2:autoregressive"),NUL<-0)
  Ans <- ifelse(option != 1, stop("The current version of the package does not supply whole sample fitting in simulation.collinear, please set 'option = 1'"), NUL <-0)

  #setting up a result matrix for lasso
  lasso.result <- matrix(0, n.resample, 5)

  #setting up a result matrix for EN
  EN.result <- matrix(0, n.resample, 5)

  #setting up a result matrix for ridge
  ridge.result <- matrix(0, n.resample, 4)


  #1 0 of the true coef vector
  true <- ifelse(coeff !=0, 1, 0)
  `%dopar%` <- foreach::`%dopar%`
  i <- NULL
  for(j in 1:n.resample){

    #resampled data
    data <- simulation.generation.data(n = n, coeff = coeff, collinear = collinear, sig = sig,
                                       split.prop = split.prop,option=option, matrix.option = matrix.option)

    a <- seq(0,1,step.alpha)

    search <- foreach::foreach(i = a, .combine = rbind) %dopar% {

      set.seed(12345);cv <- glmnet::cv.glmnet(as.matrix(data[[2]]), data[[1]], paralle = parallel, alpha = i)
      data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se],lambda.1se = cv$lambda.1se, alpha = i)
    }

    ridge.mo <- glmnet::glmnet(as.matrix(data[[2]]), data[[1]],alpha = search[1,3], lambda = search[1,2])
    pre <- predict(ridge.mo, newx=as.matrix(data[[4]]))
    pre.err <- mean((data[[3]] - pre)^2)
                          #in sample    #out sample   lambda   alpha
    ridge.result[j,] <- c(search[1,1], pre.err, search[1,2], search[1,3])


    lasso.mo <- glmnet::glmnet(as.matrix(data[[2]]), data[[1]],alpha = search[nrow(search),3], lambda = search[nrow(search),2])
    pre <- predict(lasso.mo, newx=as.matrix(data[[4]]))
    pre.err <- mean((data[[3]] - pre)^2)
    coe.lasso <- coef(lasso.mo)[-1]
    coe.lasso <- ifelse(coe.lasso != 0, 1,0)
    mo.error <- mean(abs((true-coe.lasso))^2)
                                  #in sample    #out sample    lambda               alpha
    lasso.result[j,] <- c(search[nrow(search),1], pre.err,  search[nrow(search),2],search[nrow(search),3],mo.error)


    EN <- search[c(-1,-nrow(search)),]
    EN <- EN[which(EN[,1] == min(EN[,1])),]
    EN.mo <- glmnet::glmnet(as.matrix(data[[2]]), data[[1]],alpha = EN[3], lambda = EN[2])
    pre <- predict(EN.mo, newx=as.matrix(data[[4]]))
    pre.err <- mean((data[[3]] - pre)^2)
    coe.en <- coef(EN.mo)[-1]
    coe.en <- ifelse(coe.en != 0, 1,0)
    mo.error <- mean(abs((true-coe.en))^2)
    #                     in  out     Lam   alp
    EN.result[j,] <- unlist(c(EN[1],pre.err,EN[2],EN[3],mo.error))



  }
  #Monte Carlo (Lasso)
  lasso.result.sim <- MC.result(lasso.result[,2])

  #Monte Carlo (EN)
  EN.result.sim <- MC.result(EN.result[,2])

  #Monte Carlo (ridge)
  Ridge.result.sim <- MC.result(ridge.result[,2])

  #model selection accuracy
  lasso.mod.error <- sum(lasso.result[,5]==0)/length(lasso.result[,5])
  EN.mod.error <- sum(EN.result[,5]==0)/length(EN.result[,5])

  final.result <- data.frame(c(lasso.result.sim[[1]],lasso.result.sim[[2]],lasso.mod.error,1),
                             c(EN.result.sim[[1]],EN.result.sim[[2]], EN.mod.error,mean(EN.result[,4])),
                             c(Ridge.result.sim[[1]],Ridge.result.sim[[2]], NA,0))
  rownames(final.result) <- c("Average External MSE","Standard Error","Average Model selection Accuracy","Average alpha")
  colnames(final.result) <- c("Lasso","Elastic Net","Ridge")

  itt.names <- paste("rep.",1:n.resample,sep="")

  lasso.result <- lasso.result[,1:4]
  colnames(lasso.result) <- c("(Int) Prediction Error","(Ext) Prediction Error","(2D 10-fold CV) opt.lambda","(2D 10-fold CV) opt.alpha")
  rownames(lasso.result) <- itt.names

  EN.result <- EN.result[,1:4]
  colnames(EN.result) <- c("(Int) Prediction Error","(Ext) Prediction Error","(2D 10-fold CV) opt.lambda","(2D 10-fold CV) opt.alpha")
  rownames(EN.result) <- itt.names

  colnames(ridge.result) <- c("(Int) Prediction Error","(Ext) Prediction Error","(2D 10-fold CV) opt.lambda","(2D 10-fold CV) opt.alpha")
  rownames(ridge.result) <- itt.names

  return(list("final.table"=final.result,"lasso" = lasso.result,"Elastice Net"=EN.result,"ridge" = ridge.result))
}



