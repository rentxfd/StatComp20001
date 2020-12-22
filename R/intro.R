#' @title A solution of solving equation by dichotomy in fixed interval using R
#' @description A solution of solving equation by dichotomy in fixed interval using R
#' @param f the equation to be solved
#' @param a the lower bound of interval
#' @param b the upper bound of interval
#' @param eps the required precision
#' @return the root of the equation 
#' @importFrom Rcpp evalCpp
#' @useDynLib StatComp20001
#' @examples
#' \dontrun{
#'   f <- function(x){
#'    x^3 - x - 1
#'    }
#'   Solveequation(f,1,2,1e-6)
#'   }
#' @export

Solveequation <- function(f,a,b,eps=1e-5){
  if (f(a)*f(b)>0)
    list(fail="finding root is fail!")
  else{
    repeat {
      if (abs(b-a)<eps) break      
      x <- (a+b)/2
      if (f(a)*f(x)<0)
        b <- x
      else
        a <- x
    }
    list(root=(a+b)/2, fun=f((a+b)/2))
  }
}










#' @title Find the maximum likelihood function for the normal distribution
#' @description Find the maximum likelihood function for the normal distribution
#' @param theta A parameter vector containing the mean and variance
#' @param x the sample data
#' @return maximum likelihood estimation of the mean and variance
#' @importFrom Rcpp evalCpp
#' @useDynLib StatComp20001
#' @examples
#' \dontrun{
#'   theta0 <- c(0,1) 
#'   X <- 1:50 
#'   result <- optim(theta0,norm.fun,x=X) 
#'   result
#'   }
#' @export

norm.fun <- function(theta,x) 
{	
  mu <- theta[1]
  sigma <- theta[2] 	
  n <- length(x)
  logL <- -0.5*n*log(2*pi)-0.5*n*log(sigma)-(1/(2*sigma))*sum((x-mu)**2) 
  return(-logL)  
}

  
  

#' @title Interval estimation and hypothesis testing function of mean when the variance of a normal population is known
#' @description Interval estimation and hypothesis testing function of mean when the variance of a normal population is known
#' @param x the sample data
#' @param mu mean of test
#' @param sigma known variance
#' @param a the confidence level
#' @param alternative Unilateral or bilateral inspection
#' @return interval estimation and hypothesis testing
#' @importFrom Rcpp evalCpp
#' @importFrom stats qnorm
#' @useDynLib StatComp20001
#' @examples
#' \dontrun{
#'   X <- 1:50 
#'   N.test(x=X,mu=24,sigma=9)
#'   }
#' @export

N.test<- function(x,mu,sigma,a=0.05,alternative='two.side'){  
  xbar <- mean(x)
   n <- length(x)	 

if(alternative=='two.side') 
{
  N <- qnorm(a/2) 
  area1 <- xbar-(sigma/n)**0.5*N
  area2 <- xbar+(sigma/n)**0.5*N
  area=c(area1,area2)
  NT=(xbar-mu)/(sigma/n)**0.5 
  TF='reject' 
  if ( abs(NT)<abs(N))  
  {
    TF='accept'  
  }
  result=list('accept_or_reject'=TF,'interval estimation'=area)
  
}

if(alternative=='less')
{
  N=qnorm(a)
  area1=xbar+(sigma/n)**0.5*N
  area2=Inf  
  area=c(area1,area2)
  
  NT=(xbar-mu)/(sigma/n)**0.5
  TF='reject'
  if (NT>N)
  {
    TF='accept'
  }
  result=list('accept_or_reject'=TF,'interval estimation'=area)
  
}

if(alternative=='greater')
{
  N=-qnorm(a)   
  area1=-Inf
  area2=xbar-(sigma/n)**0.5*N
  
  area=c(area1,area2)
  
  NT=(xbar-mu)/(sigma/n)**0.5
  TF='reject'
  if (NT<N)
  {
    TF='accept'
  }
  result=list('accept_or_reject'=TF,'interval estimation'=area)
}

return(result)
}




