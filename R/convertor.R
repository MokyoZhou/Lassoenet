#' Convert a cleaned dataset into a model.matrix
#'
#' This function converts a raw dataset into a user defined model.matrix.
#'
#' @param data A well-cleaned raw \code{data.frame} that you would like to convert into a model.matrix as an input for the glmnet package.
#' @param response The column location of y in the dataset.
#' @param x.indices The columns of Xs that you would like to convert into a \code{model.matrix}.
#'
#' @return a list with elements:
#'          \item{response}{the original response}
#'          \item{X.model.matrix}{the \code{model.matrix} for the predictors with the intercept column removed}
#'
#' @author Mokyo Zhou
#'
#' @examples
#'
#' data <- data.frame(matrix(rnorm(100*5),100,5))
#' #Assume y is on the first coordinate
#' #and the rest are Xs:
#' model_mat <- convertor(data = data, response = 1, x.indices = seq(2,5))
#'
#'
#' @export
#'
#'



convertor <- function(data, response, x.indices){

  Ans <- suppressWarnings(ifelse(is.na(as.numeric(response)),stop("please only input a singular value for the response coordinate"),NUL<-0))
  Ans <- ifelse((as.numeric(response) %in% seq(1,ncol(data)))==FALSE ,stop("please make sure you have the correct y coordinate"),NUL<-0)
  Ans <- ifelse(sum((as.numeric(x.indices) %in% seq(1,ncol(data)))==FALSE)!=0 ,stop("please make sure you have the correct x coordinates"),NUL<-0)
  Ans <- ifelse(sum((x.indices == floor(x.indices)) == FALSE) != 0, stop("please make sure you only enter integer X coordinates"), NUL <- 0)
  Ans <- ifelse(sum(as.numeric(x.indices) == as.numeric(response))!=0 ,stop("Please make sure Xs do not have the same coordinate as y"),NUL<-0)
  Ans <- ifelse(length(x.indices)!= length(unique(x.indices)) ,stop("Please make sure you have only entered each X coordinate once"),NUL<-0)
  xx <- data[,x.indices]
  xx <- data.frame(xx)
  colnames(xx) <- names(data)[x.indices]

  x.ma <- model.matrix(~.,xx)
  x.ma <- x.ma[,-1]

  return(list("response" = data[,as.numeric(response)], "X.model.matrix"= x.ma))

}

