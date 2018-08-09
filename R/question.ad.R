#'Internal lassoenet functions
#'
#'Internal lassoenet functions
#'
#'
#'@section Details: These are not intended for use by users.
#'This function is a part of the \code{\link{interactiver}} function.
#'This function asks a series of questions and provides pratical suggestions. It returns a collection of inputs for the Adaptive Lasso model.
#'
#'@param data A well-cleaned \code{data.frame}.
#'@param response The location of the response within the \code{data.frame}.\cr
#'@param x.indices The locations of the predictors within the \code{data.frame}.\cr
#'
#'
#'@author Mokyo Zhou
#'
#'
#'@export
#'



question.ad <- function(data = data, response = response, x.indices = x.indices){
  readline("8. The process will run the adaptive lasso for inference. The
adaptive lasso reduces bias by allowing each covariate to have a different
penalty. Thus, useful covariates can have a smaller penalty while irrelevant
coefficents will be penalised more heavily")

  cat("                                                                                        ")

  readline("9. Use either OLS or ridge coefficients to constuct a penalty
vector to allow for different penalisation. If collinearity is a problem,
then ridge regressive might be a better starting weight.The Following is
the VIFs based on the linear model:")


  #linear model
  y.linear <- data[,as.numeric(response)]
  y <- y.linear
  x.linear <- data[,x.indices]

  linear.model.data <- data.frame(y,x.linear)

  linear.model <- lm(y ~ ., data = linear.model.data)
  cat("                                                                                          ")


  a <- car::vif(linear.model)
  print(a)

  cat("                                                                                        ")

  readline("10. Any VIFs above 5 or GVIF^(1/(2*Df)) > 2 indicate collinearity.
But keep in mind that there are many more ways in which collinearity can be
checked and we are only showing 1 of such ways")

  cat("                                                                                        ")

  A11 <- readline("11. What cofficients should be used for the construction
of the penalty vector OLS or ridge? Base your answer on the outputs above
as well as any other detection teachniques and personal judgement")
  Ans <- ifelse(A11 != "OLS" & A11 != "ridge",stop("please only input OLS or ridge"),NUL<-0)
  cat(c("Type of coefficents for inital weight construction =", A11))
  cat("                                                                                         ")
  A.gam <- readline("12. Please define the gamma sequence which you would like
to try. Note this gamma sequence can have range from 0 to positive infinite.
Please enter in the form of e.g. seq(0,5) or c(0,0.5, 1,...)")
  Ans <- ifelse(is.na(A.gam),stop("Please enter a gamma sequence"),NUL <-0)
  A.gam <- eval(parse(text=A.gam))
  Ans <- ifelse(is.numeric(A.gam)==FALSE,stop("Please make sure you have numeric inputs for gamma.seq"),NUL<-0)
  Ans <- ifelse( sum(A.gam < 0)!=0 ,stop("Please make sure you have only entered positive numbers for gamma.seq"),NUL<-0)
  Ans <- ifelse(length(A.gam)!= length(unique(A.gam)) ,stop("Please make sure you have only entered each number once for gamma.seq"),NUL<-0)
  Ans <- ifelse(length(A.gam) == 1,stop("Please make sure you have more than 1 element within the gamma.seq"),NUL<-0)
  cat("The gamma sequence to use =", A.gam)

  cat("                                                                                          ")
  re.B <- readline("13. Please define the number of residual bootstrappings
that you would like to do. Please note that the more bootstrappings you do,
the less uncertain will your result be but more time consuming")
  Ans <- ifelse(is.na(re.B),stop("Please enter a single number for the number of residual bootstrappings"),NUL <-0)
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(re.B)),stop("please only input a numeric value for bootstraps"),NUL<-0))
  re.B <- as.numeric(re.B)
  Ans <- ifelse(length(re.B) != 1,stop("Please make sure you have only entered a single number for the number of residual bootstrappings"),NUL<-0)
  Ans <- ifelse(re.B < 0 ,stop("Please make sure the number of residual bootstrappings is a positive number"),NUL<-0)
  Ans <- ifelse(re.B!= floor(re.B) ,stop("Please make sure the number of residual bootstrappings is a whole number"),NUL<-0)
  Ans <- ifelse(re.B < 50 ,stop("Please make sure the number of residuals bootstrappings is bigger than 50"),NUL<-0)
  cat("Number of residual bootstrappings to do =", re.B)
  cat("                                                                                           ")

  si.B <- readline("14. Please define the significance level of the residual
bootstrap, e.g. 0.05 for a 95% confidence interval")
  Ans <- ifelse(is.na(si.B),stop("Please enter a number for the significance level"),NUL <-0)
  Ans <- suppressWarnings(ifelse(is.na(as.numeric(si.B)),stop("please only input a numeric value for significance"),NUL<-0))
  si.B <- as.numeric(si.B)
  Ans <- ifelse(length(si.B) != 1,stop("Please make sure you have only entered a single number for the significance level"),NUL<-0)
  Ans <- ifelse(si.B <= 0 || si.B >= 1,stop("Please make sure the significance level is within (0,1)"),NUL<-0)
  cat("The significance level for the residual bootstrap =", si.B)


  return(list(A11,A.gam, re.B,si.B))
}
