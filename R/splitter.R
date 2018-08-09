#' Split the data into a training set and a test set
#'
#' This function splits the \code{data.frame} into a training set and a test set with user defined proportions.
#'
#' @param data A well-cleaned \code{data.frame} that you would like to split into a training set and a test set.
#' @param splits A two-element vector with the first element specifying the training proportion and the second being the test proportion.
#' @param random If TRUE, uses a random seed so you will get a different training and testing set everytime you run this. If FALSE, the draws are based on the seed (12345678). Default is random = FALSE.
#'
#'
#' @return a list with elements:
#'          \item{train}{the training set}
#'          \item{test}{the test set}
#'
#'
#' @author Mokyo Zhou
#'
#'
#' @examples
#'
#' data <- data.frame(matrix(rnorm(100*5),100,5))
#' traintest <- datasplit(data = data, splits = c(0.8,0.2), random = FALSE)
#' training <- traintest[[1]]
#' testing <- traintest[[2]]
#'
#'
#'
#' @export
#'






datasplit <- function(data = data, splits=splits,random=FALSE){

  Ans <- ifelse(sum(splits)!=1,stop("please make sure the proportions do sum to 1"),NUL<-0)
  Ans <- ifelse(sum(splits>=0)!=2,stop("please make sure both training and testing proportions are  > 0"),NUL<-0)
  Ans <- ifelse(sum(splits<=1)!=2,stop("please make sure the proportions are less than 1"),NUL<-0)
  smp_size <- floor(splits[1] * nrow(data))
  if(random == TRUE){
    set.seed(Sys.time())
  }else{
    set.seed(12345678)
  }


  train_ind <- sample(seq(1 : nrow(data)), size=smp_size, replace=FALSE)

  train <- data[train_ind,]
  test <- data[-train_ind,]

  return(list("train"=train,"test"=test))
}





