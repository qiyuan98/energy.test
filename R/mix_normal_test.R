#' genirate a population of mix normal distribution
#'
#' this function genirate n indepently identically distributed random variables
#' @param n number of iid rv
#' @param ... mean, sd, pr of each normal population
genirate_mix_normal = function(n,...)
{
  paramatic = t(matrix(c(...),nrow=3))
  pop.N=NROW(paramatic)
  sam.N=sample(1:pop.N, replace = TRUE,size = n, prob = paramatic[,3])
  output=rnorm(n,paramatic[,1][sam.N],paramatic[,2][sam.N])
  output
}

#' computes the sample skewness coeff.
#'
#' @param x the sample to compute
sk <- function(x)
{

  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}

#' estimation the power in a mixture probability
#'
#' @param alpha the pointed value
#' @param m the number of replicates during Monte Carlo Methods for Estimation
#' @param n the size of the sample
#' @param ... mean, sd, pr of each normal population
powt=function(alpha,m,n,...)
{
  test1 <- test2 <- test3 <- numeric(m)
  #critical value for the skewness test
  cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))

  # estimate power
  for (j in 1:m) {
    x=genirate_mix_normal(n,...)

    test1[j] <- as.integer(abs(sk(x)) >= cv)
    test2[j] <- as.integer(
      shapiro.test(x)$p.value <= alpha)
    test3[j] <- as.integer(
      mvnorm.etest(x, R=200)$p.value <= alpha)
  }

  output=c(mean(test1), mean(test2), mean(test3))
}

#' estimation the power in a mixture probability
#' under alpha=0.1
#' m=300
#' n=30
#' two population
#'
#' @param mean1 the mean of the first population
#' @param sd1 the sd of the first population
#' @param mean2 the mean of the second population
#' @param sd2 the sd of the second population
pow2p.a1.2H30=function(mean1,sd1,mean2,sd2)
{
  output=matrix(numeric(44),ncol = 4)
  epsilon=seq(0,1,0.1)
  output[,1]=epsilon

  for (n in 1:11)
  {
   output[n,2:4]=powt(0.1,200,30,mean1,sd1,epsilon[n],mean2,sd2,(1-epsilon[n]))
  }

  output
}


pow2p.a1.2H30.changesd2=function(mean)
{
  M=array(numeric(19*44),dim=c(11,4,19))
  sd2v=seq(1,10,0.5)
  for (n in 1:length(sd2v))
  {
    sd2=sd2v[n]
    M[,,n]=pow2p.a1.2H30(0,1,0,sd2)
  }

  M

}

#' plot the power function of a 2 population mixed nuomal distribution
#'
#' @param A the power matric
pow2p.plot = function(A)
{
  epsilon=rep(A[,1],times=3)
  method =rep(c("skewness","S_W","energy"),each=11)
  power=as.vector(A[,2:4])
  df=data.frame(epsilon=epsilon,method=method,power=power)
  printout=ggplot(data = df, mapping = aes(x = epsilon, y = power, colour = method)) +
    geom_line()+
    ylim(0,1)+
    scale_x_continuous(breaks = seq(0,1,0.1))

  printout
}

#' plot the p.d.f of two normal distribution
#'
#' @param mean1 the mean of the first population
#' @param sd1 the sd of the first population
#' @param mean2 the mean of the second population
#' @param sd2 the sd of the second population
norm.plot=function(mean1,sd1,mean2,sd2)
{
  x=seq(-15,15,0.05)
  dn.1=dnorm(x,mean = mean1, sd = sd1, log = FALSE)
  dn.2=dnorm(x,mean = mean2, sd = sd2, log = FALSE)



  xx=rep(x,times=2)
  dn=c(dn.1,dn.2)

  ch1=paste("N(",mean1,",",sd1^2,")")
  ch2=paste("N(",mean2,",",sd2^2,")")
  pop=rep(c(ch1,ch2),each=length(x))
  df=data.frame(x=xx,dn=dn,Population=pop)

  printout=ggplot(data  = df, mapping = aes(x = xx, y = dn, colour = Population)) +
    geom_line(alpha=0.7,size=1)+
    ylim(0,1)+
    theme(axis.text.x = element_blank())

  printout
}



