# Adapted from RE2C from http://software.buhmhan.com/RE2C/index.php
# by Gabriel Hoffman
# December 19, 2019

#' Fixed effect meta-analysis for correlated test statistics
#'
#' Fixed effect meta-analysis for correlated test statistics using the Lin-Sullivan method.
#'
#' @param beta regression coefficients from each analysis
#' @param stders standard errors corresponding to betas
#' @param cor correlation matrix between of test statistics.  Default considers uncorrelated test statistics 
#'
#' @details Perform fixed effect meta-analysis for correlated test statistics using method of Lin and Sullivan (2009).  By default, correlation is set to identity matrix to for independent test statistics.  
#'
#' This method requires the correlation matrix to be symmatric positive definite (SPD).  If this condition is not satisfied, results will be NA.  If the matrix is not SPD, there is likely an issue with how it was generated. 
#'
#' However, evaluating the correlation between observations that are not pairwise complete can give correlation matricies that are not SPD.  In this case, consider running \code{Matrix::nearPD( x, corr=TRUE)} to produce the nearest SPD matrix to the input. 
#'
#' @references{
#'   \insertRef{lin2009meta}{remaCor}
#' }
#'
#' @return
#' Return values:
#' \itemize{
#' \item{\code{beta}: }{effect size}
#' \item{\code{se}: }{effect size standard error}
#' \item{\code{p}: }{p-value}
#'}
#'
#' @examples
#' library(clusterGeneration)
#' library(mvtnorm)
#' 
#' # sample size
#' n = 30
#' 
#' # number of response variables
#' m = 6
#' 
#' # Error covariance
#' Sigma = genPositiveDefMat(m)$Sigma
#' 
#' # regression parameters
#' beta = matrix(.6, 1, m)
#' 
#' # covariates
#' X = matrix(rnorm(n), ncol=1)
#' 
#' # Simulate response variables
#' Y = X %*% beta + rmvnorm(n, sigma = Sigma)
#' 
#' # Multivariate regression
#' fit = lm(Y ~ X)
#' 
#' # Correlation between residuals
#' C = cor(residuals(fit))
#' 
#' # Extract effect sizes and standard errors from model fit
#' df = lapply(coef(summary(fit)), function(a) 
#'  data.frame(beta = a["X", 1], se = a["X", 2]))
#' df = do.call(rbind, df)
#' 
#' # Run fixed effects meta-analysis, 
#' # assume identity correlation  
#' LS( df$beta, df$se)
#'  
#' # Run fixed effects meta-analysis, 
#' # account for correlation 
#' LS( df$beta, df$se, C)
#' @import stats
#' @export
LS <- function(beta, stders, cor=diag(1, length(beta)) ){

   # check arguments
   if( length(beta) != length(stders) ){
      stop("Number of test statistics and standard errors must be the same")
   }
   if( nrow(cor) != ncol(cor) ){      
      stop("Correlation matrix must be square")
   }
   if( length(beta) != nrow(cor) ){      
      stop("Number of test statistics and rows of correlation matrix must be the same")
   }
   if( ! is.numeric(beta) ){      
      stop("beta must be numeric")
   }
   if( ! is.numeric(stders) ){      
      stop("stders must be numeric")
   }
   if( any(stders <= 0) ){      
      stop("All values in stders must be positive")
   }

   ## conventional FE approach  
   V <- diag(stders) %*% cor %*% diag(stders)

   # invert V, if not valid then warn and return empty results
   # Vinv <- solve(V)
   Vinv <- tryCatch( solve(V),
    error = function(e){
      warning("Covariance matrix is not invertable. Returning NA values.")
      NA
    }
  )  
  if( length(Vinv) == 1 ){
   if( is.na(Vinv) )
      return( data.frame(beta = NA, se = NA, p = NA) )
  }

   ones <- matrix(rep(1,length(beta)),nrow=1)
   
   newx <- (ones %*% Vinv %*% beta) / (ones %*% Vinv %*% t(ones))
   newv <- 1 / (ones %*% Vinv %*% t(ones))
   
   if( newv > 0 ){    
      newstd <- sqrt(newv)  
      newp <- pchisq(newx*newx/newv, 1, lower.tail=FALSE)
   }else{
      warning("Covariance matrix is not invertable. Returning NA values.")
      newstd = NA
      newp = NA
   }

   data.frame(beta = newx, se = newstd, p = newp)
}





