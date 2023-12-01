#' The covariance of symbolic histogram variable
#'
#' @name hist_cov
#' @aliases hist_cov
#' @description This function compute the covariance of symbolic histogram variable.
#' @usage hist_cov(object, var1, var2, method)
#' @param object A MatH object.
#' @param var1 A symbolic histogram variable.
#' @param var2 A symbolic histogram variable.
#' @param method The method to calculate covariance. ex: 'BG', 'BD', 'B', 'Wass'
#' @returns Return covariance of the histogram variable.
#' @examples
#' hist_cov(Blood, 'Cholesterol', 'Hemoglobin', method = 'BG')
#' @export

hist_cov <- function(object, var1, var2, method){
  location_var1 <- which(colnames(object@M) == var1)
  location_var2 <- which(colnames(object@M) == var2)
  nr <- nrow(object@M)
  nc <- ncol(object@M)
  MatH.mean <- function(object) {
    MAT <- matrix(NA, nr, nc)
    rownames(MAT) <- rownames(object@M)
    colnames(MAT) <- colnames(object@M)
    for (i in 1:nr) {
      for (j in 1:nc) {
        if (length(object@M[i, j][[1]]@x) > 0) {
          MAT[i, j] <- object@M[i, j][[1]]@m
        }
      }
    }
    return(mat = MAT)
  }

  MatH.sd <- function(object) {
    MAT <- matrix(NA, nr, nc,
                  dimnames = list(rownames(object@M), colnames(object@M)))
    for (i in 1:nr) {
      for (j in 1:nc) {
        if (length(object@M[i, j][[1]]@x) > 0) {
          MAT[i, j] <- object@M[i, j][[1]]@s
        }
      }
    }
    return(mat = MAT)
  }

  Gj <- function(a, b, p, hmean){
    if (sum((a + b) * p) / 2 <= hmean){
      return(-1)
    } else {
      return(1)
    }
  }

  Qj <- function(a, b, hmean){
    return((a - hmean)^2 + (a - hmean) * (b - hmean) + (b - hmean)^2)
  }

  QQ <- function(a, b, hmean1, hmean2){
    return(outer((a - hmean1), (b - hmean2)))
  }

  get_pvars <- function(i){
    object1 <- object@M[i, location_var1][[1]]
    object2 <- object@M[i, location_var2][[1]]
    p1 <- object1@p[2:length(object1@p)]
    p2 <- object1@p[1:(length(object1@p) - 1)]
    p3 <- object2@p[2:length(object2@p)]
    p4 <- object2@p[1:(length(object2@p) - 1)]
    pvar1 <- p1 - p2
    pvar2 <- p3 - p4
    return(list(pvar1 = pvar1, pvar2 = pvar2))
  }

  get_GQ <- function(i){
    object1 <- object@M[i, location_var1][[1]]
    object2 <- object@M[i, location_var2][[1]]
    lenx1 <- length(object1@x)
    lenx2 <- length(object2@x)
    p <- get_pvars(i)
    Q1 <- Qj(object1@x[1:(lenx1 - 1)], object1@x[2:lenx1], hist_mean(object, var1))
    Q2 <- Qj(object2@x[1:(lenx2 - 1)], object2@x[2:lenx2], hist_mean(object, var2))
    G1 <- Gj(object1@x[1:(lenx1 - 1)], object1@x[2:lenx1], p$pvar1, hist_mean(object, var1))
    G2 <- Gj(object2@x[1:(lenx2 - 1)], object2@x[2:lenx2], p$pvar2, hist_mean(object, var2))
    return(list(Q1 = Q1, Q2 = Q2, G1 = G1, G2 = G2))
  }

  get_QQ <- function(i){
    object1 <- object@M[i, location_var1][[1]]
    object2 <- object@M[i, location_var2][[1]]
    lenx1 <- length(object1@x)
    lenx2 <- length(object2@x)
    Q1 <- QQ(object1@x[2:lenx1], object2@x[2:lenx2],
             hist_mean(object, var1), hist_mean(object, var2))
    Q2 <- QQ(object1@x[2:lenx1], object2@x[1:(lenx2 - 1)],
             hist_mean(object, var1), hist_mean(object, var2))
    Q3 <- QQ(object1@x[1:(lenx1 - 1)], object2@x[2:lenx2],
             hist_mean(object, var1), hist_mean(object, var2))
    Q4 <- QQ(object1@x[1:(lenx1 - 1)], object2@x[1:(lenx2 - 1)],
             hist_mean(object, var1), hist_mean(object, var2))
    return(list(Q1 = Q1, Q2 = Q2, Q3 = Q3, Q4 = Q4))
  }

  if (method == 'BG'){
    result <- sum(MatH.mean(object)[, location_var1] *
                    MatH.mean(object)[, location_var2])/nrow(object@M) -
      hist_mean(object, var1) * hist_mean(object, var2)
    return(result)
  } else if (method == 'BD'){
    ss <- 0
    for (i in 1:nr){
      ss <- ss + sum(get_GQ(i)$G1 * get_GQ(i)$G2 *
                       outer(get_pvars(i)$pvar1, get_pvars(i)$pvar2) *
                       outer(get_GQ(i)$Q1, get_GQ(i)$Q2)^0.5)
    }
    return(ss / (3 * nrow(object@M)))
  } else if (method == 'B'){
    ss <- 0
    for (i in 1:nr){
      ss <- ss + sum((2 * get_QQ(i)$Q1 + get_QQ(i)$Q2 + get_QQ(i)$Q3 + 2 * get_QQ(i)$Q4) *
                       outer(get_pvars(i)$pvar1, get_pvars(i)$pvar2))
    }
    return(ss / (6 * nrow(object@M)))
  } else if (method == 'Wass'){
    CM <- sum(MatH.mean(object)[, location_var1] *
                MatH.mean(object)[, location_var2])/nrow(object@M) -
      hist_mean_w(object, var1) * hist_mean_w(object, var2)

    s1 <- 0
    s2 <- 0
    for (i in 1:nr){
      s1 <- s1 + sum(rQQ(object@M[i, location_var1][[1]], object@M[i, location_var2][[1]]) *
                       MatH.sd(object)[i, location_var1] *
                       MatH.sd(object)[i, location_var2]) / nrow(object@M)
    }
    for (i in 1:nr){
      for (j in 1:nr){
        s2 <- s2 + sum(rQQ(object@M[i, location_var1][[1]], object@M[j, location_var2][[1]]) *
                         MatH.sd(object)[i, location_var1] *
                         MatH.sd(object)[j, location_var2]) / nrow(object@M)^2
      }
    }
    CV <- s1 - s2
    result <- CM + CV
    return(result)
  }
}

