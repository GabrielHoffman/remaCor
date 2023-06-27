# Gabriel Hoffman
# June 2, 2023


#' @useDynLib remaCor, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

msqrt = function(S){
  dcmp = eigen(S, symmetric=TRUE)
  # a = with(dcmp, vectors %*% diag(sqrt(values)) %*% t(vectors))
  with(dcmp, crossprod(values^.25 * t(vectors)))
}

# chol rotates axes, and causes problems!!
# Use eigen-axes
get_sample = function(C_chol, nu){

	# RSS ~ W(n,Sigma)

	# C_chol = chol(C)
	# C_i = rwishart(nu, C) / nu
	# C_i_chol = rwishartc(nu, C) / sqrt(nu)
	C_i_chol = rwishart_chol(nu, C_chol) / sqrt(nu)
	C_i_chol = msqrt(crossprod(C_i_chol))
	V_i_chol = C_i_chol / sqrt(nu)

	beta_i = V_i_chol %*% rnorm(nrow(C_chol))

	C_i_chol = rwishart_chol(nu, C_chol) / sqrt(nu)
	# C_i_chol = msqrt(crossprod(C_i_chol))
	V_i_chol = C_i_chol / sqrt(nu)

	ones = matrix(1,1,length(beta_i))
	V_i = crossprod(V_i_chol)

	# if matrix is not invertable, return NA for statistic
	if( kappa(V_i) > 1e7 ){
		return(NA)
	}

	a = solve(V_i, beta_i)
	b = solve(V_i, t(ones))

	newx <- (ones %*% a)/(ones %*% b)
	newv <- 1/(ones %*% b)
	(newx / sqrt(newv))^2
}

get_stat_samples3 = function(V, nu, n.mc.samples, seed){

	if( exists(".Random.seed") ){
		old <- .Random.seed
	  on.exit({.Random.seed <<- old})
	}
	set.seed(seed)

	ch = chol(V)
	sapply(seq(n.mc.samples), function(i) get_sample(ch, nu))
}

# Empirical p-value for a Lin-Sullivan statistic
#  
# Compute an empirical p-value for a Lin-Sullivan statistic using Monte Carlo and fitting a gamma distribution to the samples.
#  
# @param LS.stat observed Lin-Sullivan statistic reported by \code{LS()}
# @param V variance-covariance matrix
# @param nu degrees of freedom
# @param n.mc.samples number of Monte Carlo samples
#  
#' @importFrom EnvStats egamma
.LS_empirical_pvalue = function(LS.stat, V, nu = Inf, n.mc.samples=1e4, seed=1){

	stopifnot( nrow(V) == ncol(V) )
	stopifnot( all(nu > 1) )

	# sample stats assuming finite sample size
	stats = get_stat_samples3(V, nu, n.mc.samples, seed)
	stats = stats[!is.na(stats)]

	# estimate gamma parameters
	fit.g = egamma(stats)

	# compute p-value
	pgamma(LS.stat, 
		shape = fit.g$parameters[1], 
		scale=fit.g$parameters[2], lower.tail=FALSE)
}


#' Fixed effect meta-analysis for correlated test statistics
#'
#' Fixed effect meta-analysis for correlated test statistics using the Lin-Sullivan method using Monte Carlo draws from the null distribution to compute the p-value.
#'
#' @param beta regression coefficients from each analysis
#' @param stders standard errors corresponding to betas
#' @param cor correlation matrix between of test statistics.  Default considers uncorrelated test statistics 
#' @param nu degrees of freedom
#' @param n.mc.samples number of Monte Carlo samples
#' @param seed random seed so results are reproducable
#' 
#' @details The theoretical null for the Lin-Sullivan statistic for fixed effects meta-analysis is chisq when the regression coefficients are estimated from a large sample size. But for finite sample size, this null distribution is not well characterized. In this case, we are not aware of a closed from cumulative distribution function.  Instead we draw covariance matrices from a Wishart distribution, sample coefficients from a multivariate normal with this covariance, and then compute the Lin-Sullivan statistic.  A gamma distribution is then fit to these  draws from the null distribution and a p-value is computed from the cumulative distribution function of this gamma.
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
#' 	data.frame(beta = a["X", 1], se = a["X", 2]))
#' df = do.call(rbind, df)
#' 
#' # Meta-analysis assuming infinite sample size
#' # but the p-value is anti-conservative
#' LS(df$beta, df$se, C)
#' 
#' # Meta-analysis explicitly modeling the finite sample size
#' # Gives properly calibrated p-values
#' # nu is the residual degrees of freedom from the model fit
#' LS.empirical(df$beta, df$se, C, nu=n-2)
#' @export
#' @seealso \code{LS()}
LS.empirical = function(beta, stders, cor = diag(1, length(beta)), nu, n.mc.samples=1e4, seed=1){

	# compute Lin-Sullivan test-statistic
	res = LS(beta, stders, cor)

 	V = diag(stders) %*% cor %*% diag(stders)

	# Compute empirical p-value using a gamma fit to 
	# Monte Carlo samples from the null distribution
	res$p = .LS_empirical_pvalue(with(res, (beta/se)^2), V, nu, n.mc.samples, seed )

	res
}
















