#' An interactive function for building penalised models
#'
#' This is an interactive function that aims to guide users to obtain better penalised models to their own datasets through a series of questions and suggestions.
#' I try to equip each question with some useful and practical hints that would hopefully steer the user into a more appropriate modelling direction and spark research interest. The function
#' has two main emphases: prediction accuracy or making inference. For the prediction accuracy, the function mainly support the Lasso and the Elastic Net. On the other hand, the function would mainly
#' uses the Adpative Lasso model (support both "OLS" and "ridge" weightings) for inference purposes. This function is extremely fixable and has many unique features that beyond the based fitting function
#' \code{\link[glmnet]{cv.glmnet}} from the glmnet package. For more information please see the "Details" section.
#'
#' @param data A well-cleaned \code{data.frame} which will be used for modelling. The \code{data.frame} is also required to have more rows than columns.
#' @param parallel Multi-cores parallelisation is fully supported. The default is FALSE
#' @param ncores The number of cores that you would like to use in the parallel processing. This is only needed if \code{parallel = TRUE}. Alse note that the function will automatically switch off the
#' extra connections at the end of the computation.
#'
#' @section Details: Here we will briefly dicsuss the unique features that are inherent within this interactive function. For users who are more interesed in prediction accuracy, the function offers to fit both the Lasso
#' and the Elastic Net model. The function will ask a series of performance related questions and obtain the answers. These answers will then be converted into a collection of inputs that will be used within the Lasso and the Elastic Net. Example questions
#' include: Whether repeated error curves (please see functions: \code{\link{prediction_Lasso}} and \code{\link{prediction_ElasticNet}} for more details) should be used to stabilise the cross-validation process? If not, should the "lambda.min" or "lambda.1se" be used as the optimum \eqn{\lambda}?
#' Step size of the alpha grid? Whether the dataset should be into a training and a testing set? We also aimed to put some explanation or practical suggestions into each of the questions. Furthermore, in this interactive mode, the function also offers automatic perdiction performance
#' comparison between the Lasso and the Elastic Net based on one of the two metrics selected by the user: MSE or MAE. The function will then return a summary table of the best model as well as some graphical visualisation of the best model.\cr
#'
#' For users who put more emphases on making inferences,
#' the function will fit the Apdative Lasso model with inputs converted from the answers provided by the user from questions such as: whether repeated error curves are desired, what weighting method to use "OLS" or "ridge" etc. As opposed to the prediction methods, in order to achieve better inference
#' ,the function will always use the full dataset. The function will also provide the 95\% confidence intervals for the parameters from the best Adaptive Lasso model. The confidence intervals are obtained by using residual boostrapping and the \eqn{\alpha} level of the intervals is fully adjustable
#' by the user. Please run and have fun with this interactive function, and for more technical information please visit the help pages for functions \code{\link{prediction_Lasso}}, \code{\link{prediction_ElasticNet}} and \code{\link{Adlasso}}.
#'
#' @return: a list with elements:
#' \item{model.formula}{(FOR PREDICTION FOCUS ONLY). This is the formula for the best model. Note, the value of \code{lambda} and \code{alpha}(if Elastic Net is the best model) can be found in the summary table (see below) and thus, the user can reproduce the model if needed.
#' For inference focus, I omit the Adaptive Lasso model formula here. Users who are interested in
#' how to reproduce the Adaptive Lasso results, I suggest this link \url{http://ricardoscr.github.io/how-to-adaptive-lasso.html} and use the results from the sumary table (see below).\cr}
#' \item{coefficients}{(FORE PREDICTION FOCUS ONLY). This contains the point estimates from the best model. Note, coefficents from the Lasso or the Elastic Net are biased! users are advised to interpret them with caution. Partially due to this reason, confidence intervals for the estimates are not provided within the glmnet package and thus, I will laso omit them here.
#' For users who are more interested in inference, I recommond the inference path within this function which uses the Adaptive Lasso.\cr}
#' \item{summaries}{A summary table contains all the important results of the best model. Users can use the results from this table to reproduce models if desire. For more information about the results in the table, please visit \code{\link{prediction_Lasso}}, \code{\link{prediction_ElasticNet}} and \code{\link{Adlasso}}. Furthermore, for the PREDICTION results, alongside the numeric results, the function will also offers some visualisation that presents the results in a graphical way.\cr
#' 1. (non-split, non-repeat): MSE vs ln(lambda) & the solution regularisation paths of the best model\cr
#' 2. (non-split, repeat): error curves for both the Lasso and the Elastic Net & the solution regularisation paths of the best model\cr
#' 3. (split, non-repeat): MSE vs ln(lambda) on the training set & the solution regularisation paths of the best model (on train)\cr
#' 4. (split, repeat): error curves for both the Lasso and the Elastic Net (on traint) & the solution regularisation paths of the best model (on train)}
#' \item{Confidence.intervals}{(FOR INFERENCE ONLY). A table contains the Adaptive Lasso point estimates, the corresponding 100(1 - \eqn{\alpha})\% confidence intervals and the proportion of nonzero estimates within each of the boostrapped parameter vectors.\cr}
#'
#' @author Mokyo Zhou
#'
#' @examples
#' \dontrun{
#' library(glmnet)
#' data(QuickStartExample)
#' #Please NOTE: you can access "QuickStartExample" by using data.frame(y,x).
#'
#' row.samples <- sample(1:100,250,TRUE)
#' data <- data.frame(y,x)[row.samples,]
#' result <- interactiver(data = data)
#'
#' #2-cores parallel
#' result <- interactiver(data = data, parallel = TRUE, ncores = 2)
#' #NOTE the function will automatically switch off the extra connections after the computation.
#' }
#' @export
#'





interactiver <- function(data = data, parallel = FALSE, ncores = 2) {
  if(parallel == TRUE){
    cl<-parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
  }

  readline("Thank you for using our interactive process. This process
will automatically generate the desired outputs based on your input
arguments. Therefore, please read the questions carefully and ONLY
enter your inputs in the format shown. Thank you :)")

  cat("                                                                                        ")


  readline("1. please be sure to process your data appropreiately before
using this package (e.g. missing values, outliers and assumptions checking)
as dirty data would heavily influence the result quality.")

  cat("                                                                                        ")

  A6 <- readline("2. Furthermore, the dataset must be either a rectangular
data.farme or an as.matrix object with number of rows >> than the number
of columns. Are you happy to precess? Yes or No" )
  Ans <- ifelse(A6 != "Yes" & A6 != "No",stop("please only input /Yes/ or No/"),NUL<-0)

  cat(A6)

  if(A6 == "No"){
    return(a <- readline("please process the data appropriately"))

  }
  cat("                                                                                        ")
  cat("                                                                                        ")


  A0 <- readline("3. Which is the focus inference or prediction? please only
enter inference or prediction?                                           ")
  Ans <- ifelse(A0 != "prediction" & A0 != "inference",stop("please only input /prediction/ or /inference/"),NUL<-0)
  cat(A0)
  cat("                                                                                          ")
  cat("                                                                                          ")



  common.question <- questioning(data = data)




  if (A0 == "prediction"){
    prediction.pro <- penalised_pred(data = data, parallel = parallel, response = common.question[[1]],
                   x.indices = common.question[[2]], err.curves = common.question[[3]], type.lambda = common.question[[4]])
    cat("                                                                                          ")

    readline("please note that the numbers in the plot are the order of the coefficients in the model")
    if(parallel == TRUE){doParallel::stopImplicitCluster()}
    return(prediction.pro)
  }else{
    weight.method <- question.ad(data = data, response =common.question[[1]] , x.indices = common.question[[2]])
    inference.pro <- Adlasso(data = data, interactive = TRUE, weight.method = weight.method[[1]], parallel = parallel, response = common.question[[1]], x.indices = common.question[[2]],
                             err.curves = common.question[[3]], gamma.seq = weight.method[[2]], type.lambda = common.question[[4]], B.rep = weight.method[[3]],significance = weight.method[[4]])
    if(parallel == TRUE){doParallel::stopImplicitCluster()}
    return(inference.pro)
  }



}

