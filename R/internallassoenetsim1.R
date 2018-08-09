#' Internal lassoenet simulation function
#'
#' Internal lassoenet simulation function
#'
#' @param vector.result A vector of simluation results
#'
#' @section Details: This function calculates the averaged simulation prediction errors and it calls the function
#' \code{\link{SE.resample}} for standard error computation.
#'
#' @author Mokyo Zhou
#'
#' @export
#'






MC.result <- function(vector.result){

  #Average MSE
  Average <- mean(vector.result)
  #SE of average MSE
  standarderror <- SE.resample(vector.result)

  return(list(Average, standarderror[[1]] ))
}


