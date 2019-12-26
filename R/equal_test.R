#'multivariate test for equal distribution
#'
#'implement a permuation test, using energy method to calculate method for multivariate test for equal distribution
#'@param dat data should be a matrix with each row representing a point
#'@param size the sample size
#'@param group the number of the groups of the samples
#'@param dim the dimension of the point
#'@return p value from the permutation test
equal_test <- function(dat,size,group,dim,R = 199){
  total.size <- size*group
  distance <- as.matrix(dist(dat,diag = TRUE,upper = TRUE,
                             method = "minkowski",p = 1))
  d_stat_storage <- numeric(R+1)
  index <- 1:(size*group)
  d_stat_storage[1] <- d_stat(size,group,index,distance)
  for(k in 1:R){
    index <- sample(1:total.size,total.size,replace = FALSE)
    d_stat_storage[k+1] <- d_stat(size,group,index,distance)
  }
  mean(d_stat_storage > d_stat_storage[1])
}
