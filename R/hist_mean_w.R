#' The wasserstein mean of symbolic histogram variable
#'
#' @name hist_mean_w
#' @aliases hist_mean_w
#' @description This function compute the wasserstein mean of symbolic histogram variable.
#' @usage hist_mean_w(object, var)
#' @param object A MatH object.
#' @param var A symbolic histogram variable.
#' @returns Return wasserstein mean of the histogram variable.
#' @examples
#' hist_mean_w(Blood, 'Cholesterol')
#' @export

hist_mean_w <- function(object, var){
  MatH.mean <- function(object){
    nr <- nrow(object@M)
    nc <- ncol(object@M)
    MAT <- matrix(NA, nr, nc,
                  dimnames = list(rownames(object@M), colnames(object@M)))
    for (i in 1:nr) {
      for (j in 1:nc) {
        if (length(object@M[i, j][[1]]@x) > 0) {
          MAT[i, j] <- object@M[i, j][[1]]@m
        }
      }
    }
    return(mat = MAT)
  }
  location_var <- which(colnames(object@M) == var)
  Mean <- apply(MatH.mean(object), 2, mean)
  mean.df <- t(data.frame(Mean))
  row.names(mean.df) <- 'mean'
  colnames(mean.df) <- colnames(object@M)
  result <- mean.df[location_var]
  return(result)
}
