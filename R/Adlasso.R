#' Find the optimum pair (\eqn{\gamma}, \eqn{\lambda}) for an Adpative Lasso model
#'
#' This function takes in inputs defined by the user and computes the optimum \eqn{\gamma} and \eqn{\lambda} for an Adaptive Lasso model. The function is very flexible and allows for many different settings
#' such as, repeated error curves, different weighting methods and a definable \eqn{\gamma} grid. This function also fully supports multiple-cores parallelisation. The main fitting process is cv.glmnet() from the
#' package glmnet.
#'
#' @param data A well-cleaned \code{data.frame} which will be used for modelling. The \code{data.frame} is also required to have more rows than columns.
#' @param x.indices The coordinates of the predictors that you would like to model with. Please provide a vecotr of locations e.g. seq(2,6).
#' @param response The location of the response within the \code{data.frame}.
#' @param err.curves Due to the fact that the cross-validation process is random, the result can vary qutie a bit (if without a seed). In order to stabilise the CV process, the function fits a collection of
#'                   Adaptive Lasso models, each with a different \eqn{\gamma} value, multiple times (Note, only different in \eqn{\gamma}, coefficents used for building the weights are the same). Therefore, for EACH \eqn{\gamma}, we will create multiple error curves over a range of \eqn{\lambda}s and the optimum
#'                   pair (\eqn{\gamma}, \eqn{\lambda}) is the pair that has the lowest averaged error curves value (local optimum). We have an optimum pair (\eqn{\gamma}, \eqn{\lambda}) for each \eqn{\gamma} and the global optimum pair
#'                   is the pair that has the overall lowest averaged error curved value.  Note, with this setting, the process tends to be slow. Thus, we highly suggest
#'                   multiple-cores parallelisation. You can set this argument to 0 if you do not wish to stabilise the process, in which case the seed (1234567) will be used for the CV process.
#'                   A postive integer indicates the stabilisation process is desired.
#'                   For more information about how this works, please see section details below.
#' @param weight.method The method that will be used to generate the inital set of cofficients which will then be used for constructing the initial weights for the Adpative Lasso. Under the current version, two
#' methods are supplied: "OLS" or "ridge". The default is "OLS".
#' @param gamma.seq A definable \eqn{\gamma} grid with range from 0 to \eqn{\infty}. Default is \code{c(0.5,1,2)}. Considering the time consumption, this should be chosen very carefully if \code{err.curves >0} is also desired.
#' @param type.lambda Either "lambda.min" or "lambda.1se". Default is "lambda.min. Note when \code{err.curves >0}, this argument will not be used.
#' @param B.rep The number of residual bootstrappings to do for confidence intervals of the parameters. Default is 500.
#' @param significance The significance level of the confidence intervals e.g. 100(1-\eqn{\alpha})\%. Default is 0.05.\cr
#' @param interactive If you are running this function, please ALWAYS keep this argument to FALSE, which is the default.
#' @param parallel parallelisation supported,default is FALSE.
#' @return a list with elements:
#' \item{seed}{if \code{err.curves = 0}, the seed (1234567) will be used to compute \eqn{\gamma} and \eqn{\lambda}.\cr}
#' \item{number of err.curves}{if \code{err.curves >0}, this shows how many error curves there are (for each \eqn{\gamma}).\cr}
#' \item{best gamma}{if \code{err.curves = 0} (no stabilization), this is a part of the optimum pair (\eqn{\gamma}, \eqn{\lambda}) that has the lowest cross validation error (out of a single 2D grid search) with seed (1234567). This is the usual way of finding out (\eqn{\gamma}, \eqn{\lambda}).\cr
#'
#'                   if \code{err.curves} > 0 (with stabilization), this is a part of the global optimum pair (\eqn{\gamma}, \eqn{\lambda}) that has the overall lowest averaged error curves value among all the local optimum pairs.\cr}
#' \item{best lambda}{A part of the global optimum pair (\eqn{\gamma}, \eqn{\lambda}). For more information on how the \eqn{\gamma} and \eqn{\lambda} get selected, please see the detail section.\cr}
#' \item{prediction error}{if \code{err.curves = 0} and \code{weight.method = "OLS"}, (no error curves, OLS weighting method), this is the cross validation score associated with the best (\eqn{\gamma},\eqn{\lambda}) pair from the a 2D grid search on the whole data with OLS initial weights and seed (1234567).\cr
#'
#' if \code{err.curves > 0} and \code{weight.method = "OLS"}, (error curves, OLS weighting method), this is the overall lowest averaged cross-validation scores associated with the global optimum (\eqn{\gamma}, \eqn{\lambda}) pair from using the whole dataset and OLS initial weights.\cr
#'
#' if \code{err.curves = 0} and \code{weight.method = "ridge"}, (no error curves, ridge weighting method), this is the cross validation score associated with the best (\eqn{\gamma},\eqn{\lambda}) pair from the a 2D grid search on the whole data with ridge initial weights and seed (1234567).\cr
#'
#' if \code{err.curves > 0} and \code{weight.method = "ridge"}, (error curves, ridge weighting method), this is the overall lowest averaged cross-validation scores associated with the global optimum (\eqn{\gamma}, \eqn{\lambda}) pair from using the whole dataset and ridge initial weights.\cr}
#' \item{prediction_lower}{The lower bound for the prediction error. For more information see details below.\cr}
#' \item{prediction_upper}{The upper bound for the prediction error. For more information see details below.\cr}
#' \item{ridge lambda}{if \code{weight.method = "ridge"}, a 3D Cross-Validation grid search will be carried out and this is the optimal lambda for the ridge initial weights. See section details for more information.\cr}
#' \item{\%null deviance explained}{This can be seen as an indicator of goodness of fit.\cr}
#' \item{CIs}{The 100(1-\eqn{\alpha})\% confidence intervals for the parameters. The confidence intervals are constructed by using residual bootstrapping.
#' The \eqn{\alpha} level can be defined by the user. Please note that a CI of (. , .) means the algorithm failed to estimate
#' a valid CI for the corresponding coefficient. However, the proportion of non-zero estimates out of B.rep bootstraps will also be given, and thus, the user can
#' still gain some insight.\cr}
#'
#'
#' @section Details: This function further develops on the cv.glmnet() function from the glmnet package to allow for more flexibility.
#' The glmnet package itself does not directly support the fitting of Adapative Lasso models. This function wraps around the main fitting function cv.glmnet()
#' and thus, provides a direct fitting process of the Adpative Lasso model. The function, under this version of the package, offers two methods for the
#' construction of the initial weights, "OLS" and "ridge". If the "OLS" method has been selected, an \code{lm()} object will be fitted and coefficients from the fit (except for the intercept) will
#' be used to create the initial weights. If the "ridge" method has been selected, an \code{cv.glmnet(..., alpha = 0)} object will be fitted and the correpsonding coefficients (except for the intercept) will be used for building
#' the weights.
#'
#' The function also offers an alternative to compute the global optimum (\eqn{\gamma}, \eqn{\lambda}) pair by averaging across the error curves instead of using a fixed seed. More specifically, for the "OLS" method
#' , after the coefficients have been obtained from a \code{lm()} fit and are converted into the initial weights, the stabilisation process takes a double looping structure where the outer layer contains the \eqn{\gamma} grid and the inner loop builds multiple Adaptive Lasso models for each \eqn{\gamma} in the outer layer.
#' In this way, for each \eqn{\gamma}, we create say, B, Adaptive Lasso models and thus, this results in B error curves over a range of \eqn{\lambda}s. Note, the weights are the same throughout this stabilisation process. For each \eqn{\gamma} then, the
#' function finds the local optimum pair (\eqn{\gamma}, \eqn{\lambda}) by averaging across these error curves and find the pair that has the lowest averaged cross validation errors. After the function finds
#' all the local optimum pairs, the golbal optimum pair is the pair that has the overall lowest averaged cross validation error. From experneice, for medium
#' size datasets, with \code{err.curves} larger than 1500, the global optimum (\eqn{\gamma}, \eqn{\lambda}) will usually converge to stable values that consistently
#' achieves the overall lowest averaged across error curves value. This is a 2D stabilisation process.
#'
#' When the method "ridge" is selected, the inital ridge coefficents is obtained by a stabilisation process that averages arcoss the error curves (The first stabilisation). Then after the coefficients have been obtained and converted into the inital weights,
#' a 2-dimensional stabilisation process similar to the "OLS" method above will then takes place. Thus, when the "ridge" method is selected, we are stabilising a 3-dimensional process with the first dimension being the ridge coefficents recovery.
#'
#' When \code{err.curves > 0}, the 95 percent confidence interval for the prediction error (overall lowest averaged error curves value) is generated as follows:
#' from the corresponding error curves for the global optimum \eqn{\gamma}, the cross-validation scores corresponding to the global optimum (\eqn{\gamma}, \eqn{\lambda}) are extracted and the command quantile()
#' is then used to compute the 95 percent confidence interval. When we are not stabilising the process e.g.
#' \code{err.curve = 0}, we compute the (\eqn{\gamma}, \eqn{\lambda}) pair with the seed (1234567) and the associated CI is computed by using the standard error
#' provided by the glmnet package and assuming normality.
#'
#'
#' @author Mokyo Zhou
#'
#' @examples
#'
#' library(glmnet)
#' data(QuickStartExample)
#' #please NOTE: You can access "QuickStartExample" by using: data.frame(y,x).
#'
#'
#' #no error curves, weight.method="OLS", gamma.seq = c(0.5,1,2),type.lambda = "lambda.min"
#' result <- Adlasso(data = data.frame(y,x), x.indices = seq(2,21), response = 1, err.curves = 0 ,
#'                  weight.method = "OLS", gamma.seq = c(0.5,1,2), type.lambda="lambda.min")
#'
#' # 100 error curves, weight.method="OLS", gamma.seq = seq(0,5), 2-cores parallel processing
#' #cl <- parallel::makeCluster(2)
#' #doParallel::registerDoParallel(cl)
#' result <- Adlasso(data = data.frame(y,x), x.indices = seq(2,21), response = 1, err.curves = 100,
#'                  weight.method = "OLS", gamma.seq = seq(0,5), parallel = TRUE)
#'
#' # no error curves, weight.method = "ridge", gamma.seq = seq(1,10), type.lambda = "lambda.1se"
#' result <- Adlasso(data = data.frame(y,x), x.indices = seq(2,21), response = 1, err.curves = 0,
#'                  weight.method = "ridge", gamma.seq = seq(1,10), type.lambda = "lambda.1se")
#'
#' #80 error curves, weight.method = "ridge", gamma.seq=c(0,0.5,1,2,2.5),with parallel (2 cores)
#' #cl <- parallel::makeCluster(2)
#' #doParallel::registerDoParallel(cl)
#' result <- Adlasso(data = data.frame(y,x), x.indices = seq(2,21), response = 1, err.curves = 80,
#'                  weight.method = "ridge", gamma.seq = c(0,0.5,1,2,2.5), parallel = TRUE)
#'
#'
#'
#'
#' @export
#'
#'
#'
#'














Adlasso <- function(data = data, x.indices = x.indices, response = response, err.curves = 0, weight.method="OLS", gamma.seq = c(0.5,1,2),type.lambda = "lambda.min", B.rep = 500,
                    significance = 0.05, interactive =FALSE,  parallel = FALSE){


  if(dim(data)[2] > dim(data)[1]){
    stop("please make sure your dataset has more rows than columns")
  }

  Ans <- suppressWarnings(ifelse(is.na(as.numeric(response)),stop("please only input a single value for the response coordinate"),NUL<-0))
  Ans <- ifelse((as.numeric(response) %in% seq(1,ncol(data)))==FALSE ,stop("please make sure you have the correct y coordinate"),NUL<-0)
  Ans <- ifelse(length(response)!=1 ,stop("please make sure you have only entered a single location for y"),NUL<-0)
  Ans <- ifelse(is.numeric(x.indices)==FALSE,stop("Please make sure you have numeric inputs for x.indices"),NUL<-0)
  Ans <- ifelse(sum((as.numeric(x.indices) %in% seq(1,ncol(data)))==FALSE)!=0 ,stop("please make sure you have the correct x coordinates"),NUL<-0)
  Ans <- ifelse(sum(as.numeric(x.indices) == as.numeric(response))!=0 ,stop("Please make sure Xs do not have the same coordinate as y"),NUL<-0)
  Ans <- ifelse(length(x.indices)!= length(unique(x.indices)) ,stop("Please make sure you have only entered each X coordinate once"),NUL<-0)
  Ans <- ifelse(length(x.indices)== 1,stop("Please make sure you have at least two X coordinates"),NUL<-0)
  Ans <- ifelse(is.numeric(err.curves)==FALSE,stop("Please make sure you have a numeric input for err.curves"),NUL<-0)
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(err.curves)),stop("please only input a single value for the number of error curves"),NUL<-0))
  Ans <- ifelse(err.curves != floor(err.curves) ,stop("Please make sure the number of error curves is a whole number"),NUL<-0)
  Ans <- ifelse(err.curves < 0 ,stop("Please make sure the number of error curves is positive"),NUL<-0)
  Ans <- ifelse(length(err.curves)!=1 ,stop("Please make sure only enter a single number for the number of error curves"),NUL<-0)
  Ans <- ifelse(weight.method != "OLS" & weight.method != "ridge",stop("Please specify the correct weighting method"),NUL<-0)
  Ans <- ifelse(is.numeric(gamma.seq)==FALSE,stop("Please make sure you have numeric inputs for gamma.seq"),NUL<-0)
  Ans <- ifelse(sum(gamma.seq < 0)>0 ,stop("Please make sure the elements within the gamma.sep are all positive"),NUL<-0)
  Ans <- ifelse(length(gamma.seq)!= length(unique(gamma.seq)) ,stop("Please make sure you have only entered each element in the gamma.seq once"),NUL<-0)
  Ans <- ifelse(length(gamma.seq) == 1,stop("Please make sure you have more than 1 element within the gamma.seq"),NUL<-0)
  Ans <- ifelse(type.lambda != "lambda.min" & type.lambda != "lambda.1se" ,stop("Please make sure you have entered the correct lambda type"),NUL<-0)
  Ans <- ifelse(is.numeric(B.rep)==FALSE,stop("Please make sure you have a numeric input for B.rep"),NUL<-0)
  Ans <- ifelse(length(B.rep) != 1,stop("Please make sure you have only entered a single number for the number of residual bootstrappings"),NUL<-0)
  Ans <- ifelse(B.rep < 0 ,stop("Please make sure the number of residual bootstrappings is a positive number"),NUL<-0)
  Ans <- ifelse(B.rep!= floor(B.rep) ,stop("Please make sure the number of residual bootstrappings is a whole number"),NUL<-0)
  Ans <- ifelse(B.rep < 50 ,stop("Please make sure the number of residuals bootstrappings is bigger than 50"),NUL<-0)
  Ans <- ifelse(is.numeric(significance)==FALSE,stop("Please make sure you have a numeric input for significance"),NUL<-0)
  Ans <- ifelse(length(significance) != 1,stop("Please make sure you have only entered a single number for the significance level"),NUL<-0)
  Ans <- ifelse(significance <= 0 || significance >= 1,stop("Please make sure the significance level is within (0,1)"),NUL<-0)


  CI=TRUE

  data.pen <- convertor(data = data, response = response, x.indices = x.indices)

  y.linear <- data[,as.numeric(response)]
  x.linear <- data[,x.indices]

  if (weight.method == "OLS") {

    OLS.result <- OLS.approach(x.linear = x.linear, y.linear = y.linear, x = data.pen[[2]], y = data.pen[[1]],
                               parallel = parallel, dataa = data, xx.indices = x.indices, gamma.seq=gamma.seq, err.curves = err.curves, type.lambda = type.lambda, CI=CI, B.rep = B.rep,significance = significance)

    if(interactive == FALSE & err.curves ==0){
      no.OLS.adlasso <- list("seed" = OLS.result[[1]][1],"best gamma" = OLS.result[[1]][2], "best lambda" = OLS.result[[1]][3],
                             "prediction error" = OLS.result[[1]][4], "prediciton_lower" = OLS.result[[1]][5],
                             "prediction_upper" = OLS.result[[1]][6], "%null deviance explained" = OLS.result[[1]][7],
                             "CIs" = OLS.result[[2]])

      return(no.OLS.adlasso)
    }
    if(interactive ==FALSE & err.curves > 0){
      OLS.adlasso <- list("number of err.curves" = OLS.result[[1]][1],"best gamma" = OLS.result[[1]][2], "best lambda" = OLS.result[[1]][3],
                             "prediction error" = OLS.result[[1]][4], "prediciton_lower" = OLS.result[[1]][5],
                             "prediction_upper" = OLS.result[[1]][6], "%null deviance explained" = OLS.result[[1]][7],
                             "CIs" = OLS.result[[2]])

      return(OLS.adlasso)
    }



    a <- t(data.frame(OLS.result[[1]]))
    if(err.curves==0){
      rownames(a) <- "Adaptive Lasso (OLS, non-repeat)"
      colnames(a) <- c("seed","best gamma","best lambda","MSE","95% lower MSE","95% upper MSE","%null deviance explained")
    }else{
      rownames(a) <- "Adaptive Lasso (OLS, repeat)"
      colnames(a) <- c("number of error curves","best gamma","best lambda","MSE","95% lower MSE","95% upper MSE","%null deviance explained")

    }

    return(list("summaries" = a, "Confidence.intervals"=OLS.result[[2]]))

  }else if(weight.method == "ridge"){
    ridge.result <- ridge.approach(x = data.pen[[2]], y = data.pen[[1]], parallel = parallel,
                                   dataa = data, xx.indices = x.indices, gamma.seq=gamma.seq, err.curves = err.curves, type.lambda =type.lambda, CI=CI, B.rep = B.rep, significance=significance)


    if(interactive == FALSE & err.curves ==0){
      no.ridge.adlasso <- list("seed" = ridge.result[[1]][1],"best gamma" = ridge.result[[1]][2], "best lambda" = ridge.result[[1]][3],
                             "prediction error" = ridge.result[[1]][4], "prediciton_lower" = ridge.result[[1]][5],
                             "prediction_upper" = ridge.result[[1]][6], "ridge lambda" = ridge.result[[1]][7], "%null deviance explained" = ridge.result[[1]][8],
                             "CIs" = ridge.result[[2]])

      return(no.ridge.adlasso)
    }
    if(interactive ==FALSE & err.curves > 0){
      ridge.adlasso <- list("number of err.curves" = ridge.result[[1]][1],"best gamma" = ridge.result[[1]][2], "best lambda" = ridge.result[[1]][3],
                          "prediction error" = ridge.result[[1]][4], "prediciton_lower" = ridge.result[[1]][5],
                          "prediction_upper" = ridge.result[[1]][6], "ridge lambda" = ridge.result[[1]][7] , "%null deviance explained" = ridge.result[[1]][8],
                          "CIs" = ridge.result[[2]])

      return(ridge.adlasso)
    }


    b <- t(data.frame(ridge.result[[1]]))
    if(err.curves==0){
      rownames(b) <- "Adaptive Lasso (ridge, non-repeat)"
      colnames(b) <- c("seed","best gamma","best lambda","MSE","95% lower MSE","95% upper MSE", "ridge lambda", "%null deviance explained")
    }else{
      rownames(b) <- "Adaptive Lasso (ridge, repeat)"
      colnames(b) <- c("number of error curves","best gamma","best lambda","MSE","95% lower MSE","95% upper MSE",
                                       "best ridge lambda","%null deviance explained")

    }

    return(list("summaries"=b,"Confidence.intervals"=ridge.result[[2]]))
  }


}


