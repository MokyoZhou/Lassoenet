#'Simulation studies for Lasso and Adpative Lasso
#'
#'This function calls the data generating function \code{simulation.generation.data} and can be used
#'to conduct simulation studies for comparison between the Lasso and the Adaptive Lasso.
#'
#' @param n.resample The number of simulation datasets to generate\cr
#' @param n The number of rows in the each \code{data.frame}
#' @param coeff A vector a true coefficients fixed acrossed the simulated datasets\cr
#' @param matrix.option 1: Using an Exchangeable correlation matrix to simulate the predictors\cr
#'                      2: Using an Autoregressive correlation matrix to simulate the predictors\cr
#' @param collinear The correlation levels within the \code{matrix.option}
#' @param sig The model inherent error, the \eqn{\sigma^2}
#' @param option 1: split the dataset according to \code{c(split.prop, 1 - split.prop)}\cr
#'               2: Use the whole dataset. Note When \code{option = 2}, the \code{split.prop} will be ignored\cr
#'               Note: When using \code{simulation.adlasso} please make sure "option = 2" under the current package.
#' @param parallel Parallelisation
#' @return A list of elements:
#' \item{final.table}{A summary table contains: 1.The averaged number of incorrectly classified coefficients. \cr
#'                                              2.The proportion of times the each method correctly chooses the true model\cr}
#' \item{lasso}{A matrix full with simulation results for the Lasso method. Users can freely use this result matrix to obtain further insight\cr}
#' \item{Adpative Lasso}{A matrix full with simualtion results for the Adatpive Lasso method. Users can freely use this result matrix to obtain further insight\cr}
#'
#' @section Details: The function is one of the core function for the simulation studies. The function supports comparison between the Lasso and the Adaptive Lasso in terms of making inferenes.
#' This function calls the function \code{simulation.generation.data} and thus, users can study different datasets of their liking. The function uses the coefficients from the Lasso fit to construct the
#' initial weights. The function investigates and compare between three different \eqn{\gamma}s e.g. c(0.5,1,2). The function provide a summary table for the simulation results. Furthermore, the function also
#' poduces visulisation of bias, variance and mean-squared error for the sampling distributions of the coefficients. Finally, the matrix containing the simulation iterations for each of these methods is also provided. Therefore, users are free to conduct further investigation.
#'
#' @examples
#' #number of simulated dataset 20, 400 rows, vector of true coefficient is c(0.4,0.4,2.5,2.5,0,0,0,
#' #0,0,0),using exchangeable correlation matrix, correlation level is 0.2 in the autoregressive
#' #matrix,model error is 2, using the whole dataset. No parallelisation.
#' simulation1 <- simulation.adlasso(n.resample = 20, n = 400, coeff = c(0.4,0.4,2.5,2.5,rep(0,6)),
#' matrix.option = 1, collinear = 0.2, sig = 2, parallel = FALSE,option=2)
#' @author Mokyo Zhou
#'
#' @export
#'
#'



simulation.adlasso <- function(n.resample = n.resample, n = n, coeff = coeff, matrix.option = 1, collinear = collinear, sig = sig
                               ,option=2, parallel = FALSE){

  Ans <- suppressWarnings(ifelse(is.na(as.numeric(n.resample)),stop("please make sure the number of iterations is a number"),NUL<-0))
  Ans <- ifelse(length(n.resample) !=1,stop("Please make sure the number of iteratiosn is a single number"),NUL<-0)
  Ans <- ifelse(n.resample < 0,stop("please make sure the number of iterations is a positive"),NUL<-0)
  Ans <- ifelse(n.resample != floor(n.resample),stop("please make sure the number of iteration is a whole number"),NUL<-0)
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(n)),stop("please make sure the number of observation is a number"),NUL<-0))
  Ans <- ifelse(length(n) !=1,stop("Please make sure the number of observations is a single number"),NUL<-0)
  Ans <- ifelse(n < 0,stop("please make sure the number of observation is a positive"),NUL<-0)
  Ans <- ifelse(n != floor(n),stop("please make sure the number of observation is a whole number"),NUL<-0)
  Ans <- ifelse(dim(n)[1] < 5 * length(coeff),stop("please make sure you dataset has observations at least 5 times the length of your coeff"),NUL<-0)
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(coeff)),stop("please make sure all elements within coeff are numbers"),NUL<-0))
  Ans <- ifelse(length(coeff) ==1,stop("Please make sure you have more than 1 element within the coeff"),NUL<-0)
  Ans <- suppressWarnings(ifelse(sum(is.na(as.numeric(coeff)))!=0,stop("Please make sure all the elements within coeff are numbers"),NUL<-0))
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(collinear)),stop("please make sure the collinear is a number"),NUL<-0))
  Ans <- ifelse(length(collinear) != 1,stop("Please make sure you have only entered a single number for collienar"),NUL<-0)
  Ans <- ifelse(collinear < -1 || collinear > 1,stop("Please make sure collinear is within [-1,1]"),NUL<-0)
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(sig)),stop("please make sure the sig is a number"),NUL<-0))
  Ans <- ifelse(length(sig) != 1,stop("Please make sure you have only entered a single number for sig"),NUL<-0)
  Ans <- ifelse(sig < 0 ,stop("Please make sure the sig is a positive number"),NUL<-0)
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(option)),stop("please make sure the option is a number"),NUL<-0))
  Ans <- ifelse(option != 1 & option != 2,stop("Please make sure you have the correct option input, 1:split, 2: nonsplit"),NUL<-0)
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(matrix.option)),stop("please make sure the matrix.option is a number"),NUL<-0))
  Ans <- ifelse(length(matrix.option)!=1 ,stop("Please make sure you have only entered 1 number for the matrix.option input, 1:exchangeable, 2:autoregressive"),NUL<-0)
  Ans <- ifelse(matrix.option !=1 & matrix.option != 2,stop("Please make sure you have the correct matrix.option input, 1:exchangeable, 2:autoregressive"),NUL<-0)
  Ans <- ifelse(option != 2, stop("The current version of the package does not supply sample splitting in simulation.adlasso, please set 'option = 2'"), NUL <-0)
  Ans <- ifelse(length(coeff)>n, stop("please make sure you have more observations than coefficients"), NUL <-0)
  Ans <- ifelse(sum(coeff>1)<2, stop("please make sure there are at least two large effects within the coefficient vector. The algorithm does not work well when all effects are small within the coeff vector under the current version of the package. I will improve this matter in the future versions"), NUL <-0)
  Ans <- ifelse(floor(n/length(coeff))<10, stop("Please increases the sample size, otherwise the algorithm will struggle with the current conditions"), NUL <-0)


  #lasso.result
  lasso.matrix <- matrix(0,n.resample,length(coeff))

  #adaptive.lasso.result
  adlasso.matrix <- matrix(0, n.resample, length(coeff))

  #incorrect zero
  inc_zero_las <- as.numeric()
  inc_zero_ad <- as.numeric()

  #incorrect effects
  inc_effect_las <- as.numeric()
  inc_effect_ad <- as.numeric()

  #model selection accuracy
  mod.acc.las <- as.numeric()
  mod.acc.ad <- as.numeric()

  #results vector
  result.vector.las <- as.numeric()
  result.vector.ad <- as.numeric()


  #the truth
  true <-ifelse(coeff != 0, 1, 0)
  `%dopar%` <- foreach::`%dopar%`
  i <- NULL
  for (j in 1:n.resample){
    #resampling for data
    data <- simulation.generation.data(n = n, coeff = coeff, collinear = collinear, sig = sig, option = option, matrix.option = matrix.option)

    set.seed(12345); lasso.fit <- glmnet::cv.glmnet(as.matrix(data[,2:ncol(data)]),data[,1])
    linear.coefficients <- coef(lasso.fit , s=lasso.fit$lambda.1se)
    linear.coefficients <- as.vector(linear.coefficients)

    r <- c(0.5,1,2)

    search <- foreach::foreach(i = r, .combine = rbind) %dopar% {
      penalty.vector <- 1/abs(linear.coefficients[2:length(linear.coefficients)])^i


      penalty.vector[which(penalty.vector == Inf)] <- 999999999



      set.seed(12345);cv <- glmnet::cv.glmnet(as.matrix(data[,2:ncol(data)]), data[,1], parallel = parallel, penalty.factor = penalty.vector)
      data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.min],lambda.min = cv$lambda.min, gamma = i,inx = which(cv$lambda == cv$lambda.min))
    }

    #obtaining the lasso coefficients
    lasso.coef <- linear.coefficients


    #obtaining the adlasso coefficients
    #browser()
    Ad <- search
    Ad <- Ad[which(Ad[,1] == min(Ad[,1])),]
    penalty.vector <- 1/abs(linear.coefficients[2:length(linear.coefficients)])^Ad[3][[1]] #
    penalty.vector[which(penalty.vector == Inf)] <- 999999999
    Ad.mo <- glmnet::glmnet(as.matrix(data[,2:ncol(data)]), data[,1],penalty.factor = penalty.vector)
    Ad.coef <- coef(Ad.mo,s=Ad[2][[1]])
    Ad.coef <- as.matrix(Ad.coef)[,1]



    lasso.matrix[j,] <- as.vector(lasso.coef)[-1]
    lasso.ind <- ifelse(lasso.matrix[j,]!=0,1,0)

    adlasso.matrix[j,] <- as.vector(Ad.coef)[-1]
    ad.ind <- ifelse(adlasso.matrix[j,]!=0,1,0)


    inc_zero_las[j] <- sum(abs(true[which(true == 1)] - lasso.ind[which(true==1)]))
    inc_effect_las[j] <- sum(abs(true[which(true == 0)] - lasso.ind[which(true==0)]))

    inc_zero_ad[j] <- sum(abs(true[which(true == 1)] - ad.ind[which(true==1)]))
    inc_effect_ad[j] <- sum(abs(true[which(true == 0)] - ad.ind[which(true==0)]))

    mod.acc.las[j] <- ifelse(sum(abs(lasso.ind - true))!=0,1,0)
    mod.acc.ad[j] <- ifelse(sum(abs(ad.ind - true))!=0,1,0)


  }

  #plotting
  max.mat <- apply(lasso.matrix,2,max);max.mat <- max(max.mat)
  max.true <- max(coeff);Max.whole <- max(max.true,max.mat)
  min.mat <- apply(lasso.matrix,2,min);min.mat <- min(min.mat)
  min.true <- min(coeff);Min.whole <- min(min.true,min.mat)

  graphics::boxplot(lasso.matrix,main="LASSO",xlab="variable",ylab="coefficient value",ylim=c(Min.whole-0.1,Max.whole+0.1))
  lines(coeff,col=2)

  amax.mat <- apply(adlasso.matrix,2,max);amax.mat <- max(amax.mat)
  amax.true <- max(coeff);AMax.whole <- max(amax.true,amax.mat)
  amin.mat <- apply(adlasso.matrix,2,min);amin.mat <- min(amin.mat)
  amin.true <- min(coeff);AMin.whole <- min(amin.true,amin.mat)


  graphics::boxplot(adlasso.matrix,main="ADAPTIVE LASSO",xlab="variable",ylab="coefficient value",ylim=c(AMin.whole-0.1,AMax.whole+0.1))
  lines(coeff,col=2)

  #A matrix of true coefficients
  beta  <- matrix(coeff,n.resample,length(coeff),byrow=TRUE)

  #obtaining the bias
  bias1 <- colMeans(abs(lasso.matrix-beta))
  bias2 <- colMeans(abs(adlasso.matrix-beta))
  bias.whole<-c(bias1,bias2)
  bias.max <- max(bias.whole)
  bias.min <- min(bias.whole)

  #otaining the variances
  var1  <- apply(lasso.matrix,2,stats::var)
  var2  <- apply(adlasso.matrix,2,stats::var)
  var.whole<-c(var1,var2)
  var.max <- max(var.whole)
  var.min <- min(var.whole)

  #obtaining the MSE
  mse1  <- colMeans((lasso.matrix-beta)^2)
  mse2  <- colMeans((adlasso.matrix-beta)^2)
  mse.whole<-c(mse1,mse2)
  mse.max <- max(mse.whole)
  mse.min <- min(mse.whole)

  plot(bias1,type="l",xlab="Variable index",ylab="Bias",main="Bias plot",ylim=c(bias.min-0.05,bias.max+0.05))
  lines(bias2,col=2)
  graphics::mtext("Adlasso(Red) Lasso(Black)",side=3)

  plot(var2,col=2,type="l",xlab="Variable index",ylab="Variance",main="Variance plot",ylim=c(var.min-0.02,var.max+0.02))
  lines(var1,col=1)
  graphics::mtext("Adlasso(Red) Lasso(Black)",side=3)


  plot(mse1,type="l",xlab="Variable index",ylab="MSE",main="MSE plot",ylim=c(mse.min-0.03,mse.max+0.03))
  lines(mse2,col=2)
  graphics::mtext("Adlasso(Red) Lasso(Black)",side=3)

  result.data <- data.frame(c(mean(inc_zero_las + inc_effect_las),1-mean(mod.acc.las)),
                            c(mean(inc_zero_ad + inc_effect_ad),1-mean(mod.acc.ad)))

  rownames(result.data) <- c("average incorrect classifications","model accuracy")
  colnames(result.data) <- c("Lasso","Adaptive Lasso")

  X.names <- paste("x",1:dim(lasso.matrix)[2],sep="")
  itt.names <- paste("rep.",1:n.resample,sep="")
  colnames(lasso.matrix) <- X.names
  rownames(lasso.matrix) <- itt.names
  colnames(adlasso.matrix) <- X.names
  rownames(adlasso.matrix) <- itt.names


  return(list("final.table"=result.data, "lasso" = lasso.matrix, "Adaptive Lasso" = adlasso.matrix))

}




























