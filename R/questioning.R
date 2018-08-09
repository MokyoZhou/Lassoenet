#'Internal lassoenet functions
#'
#'Internal lassoenet functions
#'
#'
#'@section Details: These are not intended for use by users.
#'This function is a part of the \code{\link{interactiver}} function.
#'This function asks a series of questions and provides pratical suggestions. It returns a collection of inputs for the Lasso and the Elastic Net model.
#'
#'@param data A well-cleaned \code{data.frame}.
#'
#'@author Mokyo Zhou
#'
#'
#'@export
#'

questioning <- function (data = data){


  response <- readline("4. what is your response variable's column index.  e.g. 1                ")
  Ans <- ifelse(is.na(response),stop("Please enter a number for the response"),NUL <-0)
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(response)),stop("please only input a numeric value for the response coordinate"),NUL<-0))
  Ans <- ifelse((as.numeric(response) %in% seq(1,ncol(data)))==FALSE ,stop("please make sure you have the correct y coordinate"),NUL<-0)
  Ans <- ifelse(length(response)!=1 ,stop("please make sure you have only entered a single location for y"),NUL<-0)
  cat("response index =", response)
  cat("                                                                                        ")
  cat("                                                                                        ")

  A9 <- readline("5. what are the column indices for the covariates.  e.g.
c(seq(2,4),5,seq(6,10))                                                   ")
  Ans <- ifelse(is.na(A9),stop("Please enter the locations for the predictors"),NUL <-0)
  x.indices <- eval(parse(text=A9))
  Ans <- ifelse(is.numeric(x.indices)==FALSE,stop("Please make sure you have numeric inputs for x.indices"),NUL<-0)
  Ans <- ifelse(sum((as.numeric(x.indices) %in% seq(1,ncol(data)))==FALSE)!=0 ,stop("please make sure you have the correct x coordinates"),NUL<-0)
  Ans <- ifelse(sum(as.numeric(x.indices) == as.numeric(response))!=0 ,stop("Please make sure Xs do not have the same coordinate as y"),NUL<-0)
  Ans <- ifelse(length(x.indices)!= length(unique(x.indices)) ,stop("Please make sure you have only entered each X coordinate once"),NUL<-0)
  Ans <- ifelse(length(x.indices)== 1,stop("Please make sure you have at least two X coordinates"),NUL<-0)
  cat("covariate columns =", x.indices)
  cat("                                                                                                  ")


  #stabilise
  A2 <- readline("6. The tunning process is random. One way to stabilise
the tunning process is by averaging across multiple error curves. This might
give better solutions but more time consuming.Stabilise Yes or No?
(More info from main functions e.g. prediction_Lasso)")

  Ans <- ifelse(A2 != "Yes" & A2 != "No",stop("please only input /Yes/ or No/"),NUL<-0)
  cat("Averaging across error curves =", A2)

  if(A2 == "Yes"){
    cat("                                                                                             ")

    A3 <- readline("7. please specify the number of error curves. e.g. 100                    ")
    Ans <- suppressWarnings(ifelse(is.na(as.numeric(A3)),stop("please only input a numeric value for the error curves"),NUL<-0))
    A3 <- as.numeric(A3)
    Ans <- ifelse(length(A3)!=1,stop("Please make sure you have only input 1 single number for error curves"),NUL<-0)
    Ans <- ifelse(A3 != floor(A3) ,stop("Please make sure the number of error curves is a whole number"),NUL<-0)
    Ans <- ifelse(A3 < 0 ,stop("Please make sure the number of error curves is positive"),NUL<-0)
    cat(c("number of error curves =", A3))

    type.lambda = NULL
  }else if(A2=="No"){
    A3 <- 0

    cat("                                                                                        ")

    type.lambda <- readline("7. lambda.min or lambda.1se as the regularised parameter?
lambda.min: model has the lowest CV error but could be overfitting the data.
lambda.1se: model has CV error 1 standard error away from the lowest,usually
preferred but might not predict as good")
    Ans <- ifelse(type.lambda != "lambda.min" & type.lambda != "lambda.1se",stop("please only input /lambda.min/ or /lambda.1se/"),NUL<-0)
    cat(ifelse(type.lambda=="lambda.min",paste("Using lambda.min"),paste("Using lambda.1se")))



  }


  cat("                                                                                                    ")



  return(list(response, x.indices, A3,type.lambda))
}

