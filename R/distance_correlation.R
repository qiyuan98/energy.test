#' compute the distance covariance of two samples
#'
#' this function generate the statistic of the test for dependence based on distance
#' @param x one of the sample ,can be a vector,a matrix
#' @param y the other sample,the dimension must agree with the x
dCov <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  m <- nrow(y)
  if (n != m || n < 2)
    stop("Sample sizes must agree")
  if (! (all(is.finite(c(x, y)))))
    stop("Data contains missing or infinite values")
  Akl <- function(x) {
    d <- as.matrix(dist(x))
    m <- rowMeans(d)
    M <- mean(d)
    a <- sweep(d, 1, m)
    b <- sweep(a, 2, m)
    return(b + M)
  }
  A <- Akl(x)
  B <- Akl(y)
  dCov <- sqrt(mean(A * B))
  dCov
}
#' compute the distanc correlation of two samples which known as Rn
#'
#' this function based on the dcov and give the final result
#' @param z the matrix of the two samples
#' @param ix the index of one of sample in the matrx
#' @param dims the dimension of the matrix z
ndCov2 <- function(z, ix, dims) {
  #dims contains dimensions of x and y
  p <- dims[1]
  q1 <- dims[2] + 1
  d <- p + dims[2]
  x <- z[ , 1:p] #leave x as is
  y <- z[ix, q1:d] #permute rows of y
  return(nrow(z) * dCov(x, y)^2)
}

