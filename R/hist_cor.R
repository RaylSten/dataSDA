#' The correlation of symbolic histogram variable
#'
#' @name hist_cor
#' @aliases hist_cor
#' @description This function compute the correlation of symbolic histogram variable.
#' @usage hist_cor(object, var1, var2, method)
#' @param object A MatH object.
#' @param var1 A symbolic histogram variable.
#' @param var2 A symbolic histogram variable.
#' @param method The method to calculate correlation ex: 'BG', 'BD', 'B', 'Wass'
#' @returns Return covariance of the histogram variable.
#' @examples
#' hist_cor(Blood, 'Cholesterol', 'Hemoglobin', method = 'BG')
#' @export

hist_cor <- function(object, var1, var2, method){
  if (method == 'Wass'){
    result <- hist_cov(object, var1, var2, method = method) /
      (sqrt(hist_var_w(object, var1)) * sqrt(hist_var_w(object, var2)))
  } else{
    result <- hist_cov(object, var1, var2, method = method) /
      (sqrt(hist_var(object, var1)) * sqrt(hist_var(object, var2)))
  }
  return(result)
}



