#' Forest plot of coefficients 
#'
#' Forest plot of coefficients 
#'
#' @param beta regression coefficients from each analysis
#' @param stders standard errors corresponding to betas 
#'
#' @return Forest plot of effect sizes and standard errors
#'
#' @examples
#' # Generate effects
#' library(mvtnorm)
#' library(clusterGeneration )
#' 
#' n = 4
#' Sigma = cov2cor(genPositiveDefMat(n)$Sigma)
#' beta = t(rmvnorm(1, rep(0, n), Sigma))
#' stders = rep(.1, n)	
#'
#' # set names
#' rownames(Sigma) = colnames(Sigma) = LETTERS[1:n]
#' rownames(beta) = names(stders) = LETTERS[1:n]
#' 
#' # Run random effects meta-analysis,
#' # account for correlation 
#' RE2C( beta, stders, Sigma)
#'
#' # Make plot
#' plotForest( beta, stders )
#'
#' @import stats
#' @import ggplot2
#' @export
plotForest = function( beta, stders){

	# check arguments
	if( length(beta) != length(stders) ){
		stop("Number of test statistics and standard errors must be the same")
	}
	# if( nrow(cor) != ncol(cor) ){      
	# 	stop("Correlation matrix must be square")
	# }
	# if( length(beta) != nrow(cor) ){      
	# 	stop("Number of test statistics and rows of correlation matrix must be the same")
	# }
	if( ! is.numeric(beta) ){      
		stop("beta must be numeric")
	}
	if( ! is.numeric(stders) ){      
		stop("stders must be numeric")
	}

	# extract names
	if( is(beta, "matrix") ){
		ID = rownames(beta)
	}else if( is(beta, "vector") ){
		ID = names(beta)
	}
	if( is.null(ID) ){
		stop("Cannot extract rownames() or names() from beta")
	}
	# if( ! identical(ID, rownames(cor)) ){
	# 	stop("names of beta doesn't match rownames(cor)")
	# }
	# if( ! identical(ID, colnames(cor)) ){
	# 	stop("names of beta doesn't match colnames(cor)")
	# }

	df = data.frame(ID, beta, stders)
	df$ID = factor(df$ID, rev(ID))

	ggplot(df, aes(ID, beta)) + 
		theme_bw(15) + 
		geom_hline(yintercept=0, linetype="dashed", color="grey50", size=1.2) + 
		ylab(bquote(beta)) + 
		xlab("") + 
		geom_errorbar(aes(ymin=beta-1.96*stders, ymax=beta+1.96*stders), width=.1) + 
		geom_point(color="dodgerblue", size=3) + 
		coord_flip() + 
		theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5))
}

#' Correlation plot
#'
#' Correlation plot
#'
#' @param cor correlation matrix between of test statistics.  Default considers uncorrelated test statistics   
#'
#' @return Plot of correlation matrix
#'
#' @examples
#' # Generate effects
#' library(mvtnorm)
#' library(clusterGeneration )
#' 
#' n = 4
#' Sigma = cov2cor(genPositiveDefMat(n)$Sigma)
#' beta = t(rmvnorm(1, rep(0, n), Sigma))
#' stders = rep(.1, n)	
#'
#' # set names
#' rownames(Sigma) = colnames(Sigma) = LETTERS[1:n]
#' rownames(beta) = names(stders) = LETTERS[1:n]
#' 
#' # Run random effects meta-analysis,
#' # account for correlation 
#' RE2C( beta, stders, Sigma)
#'
#' # Make plot
#' plotCor( Sigma )
#'
#' @import stats
#' @import ggplot2
#' @import grid
#' @importFrom reshape2 melt
#' @export
plotCor = function( cor ){

	Var1 = Var2 = value = NA # Pass R CMD check

	ID = rownames(cor)

	df2 = melt(cor)
	df2$Var1 = factor(df2$Var1, ID)
	df2$Var2 = factor(df2$Var2, rev(ID))

	ggplot(df2, aes(Var1, Var2, fill=value)) + 
		geom_tile() + 
		scale_fill_gradient2(name = "Correlation", low="blue", mid="white", high="red", limits=c(-1, 1)) +
		xlab("") + 
		ylab("") + 
		theme_minimal() + 
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1, axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5)) 
}




