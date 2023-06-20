# Gabriel Hoffman
# June 13, 2023

#' Hottelling T^2 test for multivariate regression
#' 
#' Hottelling T^2 test compares estimated regression coefficients to specified values under the null.  This tests a global hypothesis for all specified coefficients.  It uses the F-distribution as the null for the test statistic which is exact under finite sample size.    
#' 
#' @param beta regressioin coefficients
#' @param Sigma covariance matrix of regression coefficients
#' @param n sample size used for estimation
#' @param mu_null the values of the regression coefficients under the null hypothesis.  Defaults to all zeros
#' 
#' @details The Hotelling T2 test is not defined when n - p < 1.  Returns \code{data.frame} with \code{stat = pvalue = NA}.

#' @examples
#' library(clusterGeneration)
#' library(mvtnorm)
#' 
#' # sample size
#' n = 30
#' 
#' # number of response variables
#' m = 2
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
#' # extract coefficients and covariance
#' # corresponding to the x variable
#' beta = coef(fit)['X',]
#' S = vcov(fit)[c(2,4), c(2,4)]
#' 
#' # perform Hotelling test
#' hotelling(beta, S, n)
#' 
#' @export
hotelling = function(beta, Sigma, n, mu_null = rep(0, length(beta))){

	beta = as.vector(beta)
	stopifnot(nrow(Sigma) == ncol(Sigma))
	stopifnot(nrow(Sigma) == length(beta))

	p = length(beta)

	if( n - p < 1){
		stat = pvalue = NA
		warning("Hotelling test not defined when n - p < 1")
	}else{
		# Hotelling T^2 statistic
		tsq = crossprod(beta - mu_null, solve(Sigma, beta - mu_null))

		# compute test statistic
		stat = (n-p) / (p*(n-1)) * tsq[1]

		# Compare test statistic to finite-sample null
		pvalue = pf(stat, p, n-p, lower.tail=FALSE)
	}

	data.frame(n = n, 
				p = p,
				tsq = tsq,
				stat = stat,
				p.value = pvalue)
}

