d_stat <- function(ssize = ssize,g = group,index,distance){
  total.s <- ssize * g
  ind <- index;ind <- matrix(ind,byrow = TRUE,nrow = group)
  add.w <- numeric(g)
  for (i in 1:g) {
    add.w[i] <- ssize/2 * mean(distance[ind[i,],ind[i,]])
  }
  w <- sum(add.w)
  s <- 0
  for(i in 1:(g-1)){
    for(j in (i+1):g){
      s = s + ssize^2/total.size *
        ( 2*mean(distance[ind[i,],ind[j,]]) -
            mean(distance[ind[j,],ind[j,]]) -
            mean(distance[ind[i,],ind[i,]]) )
    }
  }
  s/(g - 1)/(w/(total.s - g))
}
