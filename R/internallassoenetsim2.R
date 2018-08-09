#' Internal lassoenet simulation function
#'
#' Internal lassoenet simulation function
#'
#' @param vector.result A vector of simluation results
#'
#' @section Details: This function calculates the standard error.
#'
#' @author Mokyo Zhou
#'
#' @export
#'



SE.resample <- function(vector.result){
  result <- as.numeric()

  for(i in 1:500){
    result[i] <- mean(sample(vector.result,length(vector.result),T))
  }


  return(list(stats::sd(result)))
}
