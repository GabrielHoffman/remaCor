# Gabriel Hoffman
# June 2, 2023

# Monte Carlo simulation of Lin-Sullivan statistics given covariance and degrees of freedom
#' @importFrom mvtnorm rmvnorm
get_stat_samples = function(V, nu = rep(Inf, nrow(V)), n.mc.samples=1e5, seed=1 ){

	# use specified seed internally, then reset to original value
	old <- .Random.seed
    on.exit({.Random.seed <<- old})
    set.seed(seed)

	p = nrow(V)

	# Simulate beta from null distribution
	# multivariate normal
	beta_null = rmvnorm(n.mc.samples, rep(0, p), V)

	# Simulate from Student t
	# a ~ N(0, Sigma)
	# u ~ chisq(nu)
	# x = sqrt(nu/u) * a is MVT
	s = lapply( seq(p), function(i){
		# normal
		if( is.infinite(nu[i]) ) rep(1, n.mc.samples)
		# student t scaling
		else sqrt(nu[i] / rchisq(n.mc.samples, nu[i]))
		})
	s = do.call(cbind, s)

	# scale beta_null to have Student-t distribution
	beta_null = s * beta_null

	# Compute Lin-Sullivan statistics
	ones <- matrix(rep(1,p),nrow=1)

	Vinv.ones = solve(V, t(ones))

	numerator = as.numeric(beta_null %*% Vinv.ones)
	denom = crossprod(Vinv.ones, t(ones))[1]

	newx = numerator / denom
	newv = 1/denom
	stat = newx * newx/newv
	
	stat
}


# Empirical p-value for a Lin-Sullivan statistic
#  
# Compute an empirical p-value for a Lin-Sullivan statistic using Monte Carlo and fitting a gamma distribution to the samples.
#  
# @param LS.stat observed Lin-Sullivan statistic reported by \code{LS()}
# @param V variance-covariance matrix
# @param nu array of degrees of freedom values, one for each coefficient
# @param n.mc.samples number of Monte Carlo samples
#  
# @importFrom MASS fitdistr
#' @importFrom EnvStats egamma
.LS_empirical_pvalue = function(LS.stat, V, nu = rep(Inf, nrow(V)), n.mc.samples=1e5, seed=1){

	stopifnot( nrow(V) == ncol(V) )
	stopifnot( nrow(V) == length(nu) )
	stopifnot( all(nu > 1) )

	# sample stats assuming finite sample size
	stats = get_stat_samples(V, nu, n.mc.samples, seed)

	# # slower estimation of gamma paramters
	# fit.g = suppressWarnings(fitdistr(stats, "gamma"))
	# pgamma(LS.stat, shape = fit.g$estimate[1], rate=fit.g$estimate[2], lower.tail=FALSE)

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
#' @param nu array of degrees of freedom values, one for each coefficient
#' @param n.mc.samples number of Monte Carlo samples
#' @param seed random seed so results are reproducable
#' 
#' @details The theorical null for the Lin-Sullivan statistic for fixed effects meta-analysis is chisq when the regression coefficients being tested are normally distributed. This works  asymptotically when the regression coefficients are estimated from a sufficiently large sample.  But when the coefficients are estimated from a small sample, they have a Student-t null distrubtion that is not well approximated by a normal.  In this case, no cumulative distribution function exists in closed form.  Here we simulate coefficients from the student-t under the null given the covariance and degrees of freedom, compute the Lin-Sullivan statistic and then fit a gamma distribution to these null samples.  The p-value for the observed Lin-Sullivan statistic is then computed using this gamma approximation.
#'
#' @examples
#' library(clusterGeneration)
#' library(mvtnorm)
#' 
#' p = 6
#' C = cov2cor(genPositiveDefMat(p)$Sigma)
#' stders = rep(2, p)
#' V <- diag(stders) %*% C %*% diag(stders)
#' beta = t(rmvnorm(1, rep(5, p), V))
#' 
#' # Run fixed effects meta-analysis, 
#' # and get theoretical p-value  
#' res = LS( beta, stders, C)
#' res
#' 
#' # Get p-value from gamma fit to Monte Carlo samples 
#' # from null distrubtion
#' # This matches the theoretical p-value when betas 
#' # are drawn from a normal distribution
#' LS.empirical(beta, stders, C)
#' 
#' # Compute empirical p-value assuming 
#' # varying and finite degrees of freedom
#' nu = seq(5, p+4)
#' LS.empirical(beta, stders, C, nu)
#' @export
#' @seealso \code{LS()}
LS.empirical = function(beta, stders, cor = diag(1, length(beta)), nu = rep(Inf, length(beta)), n.mc.samples=1e5, seed=1){

	# compute Lin-Sullivan test-statistic
	res = LS(beta, stders, cor)

 	V = diag(stders) %*% cor %*% diag(stders)

	# Compute empirical p-value using a gamma fit to 
	# Monte Carlo samples from the null distribution
	res$p = .LS_empirical_pvalue(with(res, (beta/se)^2), V, nu, n.mc.samples, seed )

	res
}
















