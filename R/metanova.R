# by Gabriel Hoffman
# January 20, 2020

#' Class metanova
#'
#' Class \code{metanova} 
#'
#' @name metanova-class
#' @rdname metanova-class
#' @exportClass metanova
setClass("metanova", contains="data.frame", slots=c( method='character'))

#' show for metanova
#'
#' show for metanova
#'
#' @param object metanova object
#' @param ... other arguments, currently ignored
#'
#' @return results of show
#' @export
#' @rdname show-method
#' @aliases show,metanova-method
setMethod("show", "metanova",
function( object){
	print( object )
})

#' print for metanova
#'
#' print for metanova
#'
#' @param x metanova object
#' @param ... other arguments, currently ignored
#'
#' @return results of print
#' @export
#' @rdname print-method
#' @aliases print,metanova-method
setMethod("print", "metanova",
function( x, ...){

	txt = switch( x@method, 
		"LS" = '\n\tLin-Sullivan (fixed) meta-analysis\n\n',
		"RE2C" = '\n\tRE2C (random) meta-analysis\n\n')

	cat(txt)
	print( data.frame(x) )
})



#' Meta-analysis on multivariate regression
#'
#' Meta-analysis on multivariate regression using either LS (fixed) or RE2C (random test)
#'
#' @param fit multivariate regression object from lm()
#' @param method specify the test: "LS" or "RE2C"
#' @examples
#' # simulate data
#' N = 100
#' Y = cbind(rnorm(N), rnorm(N))
#' x = rnorm(N)
#' 
#' # fit multivariate regression
#' fit = lm(Y ~ x)
#' 
#' # evaluate LS test
#' metanova(fit, "LS")
#' 
#' # evaluate RE2C test
#' metanova(fit, "RE2C")
#' 
#' # standard anova test
#' anova(fit)
#' 
#' @import Rdpack
#' @export
metanova = function( fit, method=c("LS", "RE2C") ){

	method = match.arg( method )

	# perform data
	res = switch( method, 
		"LS" = .htls( fit ),
		"RE2C" = .htre2c( fit )[,seq(1,3)])

	# create object
	new("metanova", res, method=method)
}

# b = factor(sample(1:3, N, replace=TRUE))
# Y = cbind(rnorm(N)+as.numeric(b)/30, rnorm(N)+as.numeric(b)/30)
# fit = lm(Y ~ b)
# anova(fit)

# sumlog(metanova(fit)$p[-1])



# Hypothesis test for multivariante linear regression
# Using LS test
.htls = function( fit ){

	# for each coefficient
	res = lapply( rownames(coef(fit)), function(coef){
		# get variance-covariance matrix
		V = vcov(fit)

		# get indecies corresponding to specified coefficient
		idx = colnames(V) == paste0(':',coef)

		# get standard error from vcov matrix
		stderr = sqrt(diag(V[idx,idx]))

		# convert vcov matrix to correlation
		C = cov2cor(V[idx,idx])

		# perform meta-analysis
		# This identical to: 
		# LS( coef(fit)[coef,], c(1,1), V[idx,idx] )
		LS( coef(fit)[coef,], stderr, C )
	})
	res = do.call("rbind", res)
	rownames(res) = rownames(coef(fit))

	res
}


# Hypothesis test for multivariante linear regression
# Using RE2C test
.htre2c = function( fit ){

	# for each coefficient
	res = lapply( rownames(coef(fit)), function(coef){
		# get variance-covariance matrix
		V = vcov(fit)

		# get indecies corresponding to specified coefficient
		idx = colnames(V) == paste0(':',coef)

		# get standard error from vcov matrix
		stderr = sqrt(diag(V[idx,idx]))

		# convert vcov matrix to correlation
		C = cov2cor(V[idx,idx])

		# perform meta-analysis
		# This identical to: 
		# LS( coef(fit)[coef,], c(1,1), V[idx,idx] )
		RE2C( coef(fit)[coef,], stderr, C )
	})
	res = do.call("rbind", res)
	rownames(res) = rownames(coef(fit))

	colnames(res)[3] = 'p'

	res
}


