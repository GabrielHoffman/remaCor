# Adapted from RE2C from http://software.buhmhan.com/RE2C/index.php
# by Gabriel Hoffman
# December 19, 2019

# Log-likelihood function
#
# Log-likelihood function
#
# @param x x
# @param beta beta
# @param sigma sigma
# 
#' @importFrom mvtnorm dmvnorm
ll.fun = function(x,beta,sigma){
  -sum(dmvnorm(beta, mean=rep(x[1],length(beta)),sigma=sigma+diag(x[2],length(beta)),log=TRUE))
}

LL.fun = function(x,n,xis,etasqs,lambdas) {
  0.5*(n*log(2*pi)+sum(log(xis+x))+sum(etasqs/(lambdas+x)))
}

cal_FEprobs = function(FEs){
   pchisqs<-pchisq(FEs,df=1,ncp=0,lower.tail=FALSE)
   inst=1
   FEprobs=NULL
   for (i in seq(1, length(FEs)) ){
      FEprobs <- c(FEprobs,inst-pchisqs[i])
      inst=pchisqs[i]
   }
   return(FEprobs)
}

# Define function
#
# Define function
#
# @param nstudy number of datasets
#
#' @importFrom compiler cmpfun
#' @import stats methods
test_n = function(nstudy){
  if (nstudy == 2| nstudy == 3){
    NR_R <- function(x,n,xis,etasqs,lambdas){
      LL.fun <- function(x) {
        0.5*(n*log(2*pi)+sum(log(xis+rep(x,n)))+sum(etasqs/(lambdas+rep(x,n-1))))
      }
      optim.rst <- optim(par=c(x), fn=LL.fun, method = "L-BFGS-B",
                         lower = 0, upper=Inf)
      return(optim.rst$par[1])
    }
  } else {
      NR_R <- function(x,n,xis,etasqs,lambdas){
         init <- x
         i = 1
         while(abs(converge <- ( -0.5 * sum( 1/(xis + rep(x,n)) ) + 0.5 * sum (  etasqs / (lambdas + rep(x,n-1))^2 ) )) > 10^-8 )   {
            newx <- x - (converge / ( 0.5 * sum(  1 / (xis + rep(x,n))^2) - sum ( etasqs / (lambdas + rep(x,n-1))^3 ))  )
            x <- newx
            i = i+1
            if(i%%100==0) {
              # warning("NR failed to converge. Use optim to estimate tau^2")
              LL.fun <- function(x) {
                0.5*(n*log(2*pi)+sum(log(xis+rep(x,n)))+sum(etasqs/(lambdas+rep(x,n-1))))
              }
              optim.rst <- optim(par=c(x), fn=LL.fun, method = "L-BFGS-B",
                                 lower = 0, upper=Inf)
              return(optim.rst$par[1])
              break
            }
         }
         if(abs(init-x) > 10^4) {
            return(init)
         } else if (x < 0) {
            return(0)
         } else {
            return(x)
         }
      }
   }
   # cmpfun(NR_R,options = list(optimize=2,suppressAll=T))
   NR_R
}

