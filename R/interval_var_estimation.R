#' The estimation of symbolic variance for interval variable
#'
#' @name interval_var_estimation
#' @aliases interval_var_estimation
#' @description This function compute the estimation of symbolic variance for interval variable.
#' @usage interval_var_estimation(x, n, alpha)
#' @param x A symbolic interval variable.
#' @param n The number of elements in each group.
#' @param alpha The α-correction term related to the asymptotic mean of order statistics.
#' @returns Return mean of the interval variable.
#' @importFrom stats qnorm
#' @examples
#' data("Mushroom")
#' interval_var_estimation(Mushroom$Pileus.Cap.Width, n = 3, alpha = 3/8)
#' @export


interval_var_estimation <- function(x, n, alpha = 3/8){
  if (mode(x) == 'list'){
    numerator <- (Im(sum(apply(x, 1, sum))) - Re(sum(apply(x, 1, sum)))) / nrow(x)
  } else if (mode(x) == 'complex'){
    numerator <- (sum(data.frame(x)$max - data.frame(x)$min)) / nrow(data.frame(x))
  } else {
    print('Please enter interval variable.')
  }
  denominator <- stats::qnorm((n - alpha) / (n - 2 * alpha + 1)) - stats::qnorm((1 - alpha) / (n - 2 * alpha + 1))
  variance <- (numerator / denominator) ^ 2
  return(variance)
}
