#' The mean of symbolic histogram variable
#'
#' @name hist_mean
#' @aliases hist_mean
#' @description This function compute the mean of symbolic histogram variable.
#' @usage hist_mean(object, var)
#' @param object A MatH object.
#' @param var A symbolic histogram variable.
#' @returns Return mean of the histogram variable.
#' @examples
#' hist_mean(Blood, 'Cholesterol')
#' @export

hist_mean <- function(object, var){
  location_var <- which(colnames(object@M) == var)
  nr <- nrow(object@M)
  nc <- ncol(object@M)
  Sample_mean <- c()
  for (i in 1:nr){
    for (j in 1:nc){
      p1 <- object@M[i, j][[1]]@p[2:length(object@M[i, j][[1]]@p)]
      p2 <- object@M[i, j][[1]]@p[1:length(object@M[i, j][[1]]@p) - 1]
      p <- p1 - p2
      m2 <- c()
      for (k in 1:length(p)){
        m1 <- (object@M[i, j][[1]]@x[k] + object@M[i, j][[1]]@x[k + 1])*p[k]
        m2 <- c(m2, m1)
      }
      m3 <- sum(m2)
      Sample_mean <- c(Sample_mean, m3)
    }
  }
  Mean <- matrix(Sample_mean, nrow = nr, byrow = T,
                 dimnames = list(rownames(object@M), colnames(object@M)))
  Sample.means <- apply(Mean, 2, sum)/(2 * nr)
  mean.df <- t(data.frame(Sample.means))
  row.names(mean.df) <- 'mean'
  colnames(mean.df) <- colnames(object@M)
  result <- mean.df[location_var]
  return(result)
}
