#' The estimation of symbolic mean for interval variable
#'
#' @name interval_mean_estimation
#' @aliases interval_mean_estimation
#' @description This function compute the estimation of symbolic mean for interval variable.
#' @usage interval_mean_estimation(x)
#' @param x A symbolic interval variable.
#' @returns Return mean of the interval variable.
#' @examples
#' data("Mushroom")
#' interval_mean_estimation(Mushroom$Pileus.Cap.Width)
#' @export

interval_mean_estimation <- function(x){
  if (mode(x) == 'list'){
    return((Im(sum(apply(x, 1, sum))) + Re(sum(apply(x, 1, sum)))) / (2 * nrow(x)))
  } else if (mode(x) == 'complex'){
    return(sum(apply(data.frame(x), 1, sum)) / (2 * nrow(data.frame(x))))
  } else {
    print('Please enter interval variable.')
  }
}
