# Adapted from RE2C from http://software.buhmhan.com/RE2C/index.php
# by Gabriel Hoffman
# December 19, 2019

#' Local environment
#' 
#' Local environment
# @export
pkg.env <- new.env()


#' Random effect meta-analysis for correlated test statistics
#'
#' Random effect meta-analysis for correlated test statistics using RE2C
#'
#' @param beta regression coefficients from each analysis
#' @param stders standard errors corresponding to betas
#' @param cor correlation matrix between of test statistics.  Default considers uncorrelated test statistics 
#' @param twoStep Apply two step version of RE2C that is designed to be applied only after the fixed effect model.
#'
#' @details Perform random effect meta-analysis for correlated test statistics using RE2 method of Han and Eskin (2011), or RE2 for correlated test statistics from Han, et al., (2016).  Also uses RE2C method of Lee, Eskin and Han (2017) to further test for heterogenity in effect size. By default, correlation is set to identity matrix to for independent test statistics.
#'
#' This method requires the correlation matrix to be symmatric positive definite (SPD).  If this condition is not satisfied, results will be NA.  If the matrix is not SPD, there is likely an issue with how it was generated. 
#'
#' However, evaluating the correlation between observations that are not pairwise complete can give correlation matricies that are not SPD.  In this case, consider running \code{Matrix::nearPD( x, corr=TRUE)} to produce the nearest SPD matrix to the input. 
#'
#' @references{
#'   \insertRef{lee2017increasing}{remaCor}
#' 
#'   \insertRef{han2016general}{remaCor}
#' 
#'   \insertRef{han2011random}{remaCor}
#' }
#'
#' @return
#' Return values:
#' \itemize{
#' \item{\code{stat1}: }{statistic testing effect mean}
#' \item{\code{stat2}: }{statistic testing effect heterogeneity}
#' \item{\code{RE2Cp}: }{RE2 p-value accounting for correlelation between tests}
#' \item{\code{RE2Cp.twoStep}: }{two step RE2C test after fixed effect test.  Only evaluated if twoStep==TRUE}
#' \item{\code{QE}: }{test statistic for the test of (residual) heterogeneity}
#' \item{\code{QEp}: }{p-value for the test of (residual) heterogeneity}
#' \item{\code{Isq}: }I^2 statistic
#' }
#' \code{QE}, \code{QEp} and \code{ISq} are only evaluted if correlation is diagonal
#'
#' @examples
#' # Generate effects
#' library(mvtnorm)
#' library(clusterGeneration )
#' 
#' set.seed(1)
#' n = 4
#' Sigma = cov2cor(genPositiveDefMat(n)$Sigma)
#' beta = t(rmvnorm(1, rep(0, n), Sigma))
#' stders = rep(.1, n)
#' 
#' # Run fixed effects meta-analysis, 
#' # assume identity correlation  
#' LS( beta, stders)
#' 
#' # Run random effects meta-analysis,
#' # assume identity correlation  
#' RE2C( beta, stders)
#' 
#' # Run fixed effects meta-analysis, 
#' # account for correlation 
#' LS( beta, stders, Sigma)
#' 
#' # Run random effects meta-analysis,
#' # account for correlation 
#' RE2C( beta, stders, Sigma)
#'
#' @import stats RUnit clusterGeneration Rdpack
#' @export
RE2C <- function(beta, stders, cor=diag(1,length(beta)), twoStep = FALSE) {

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

  # Added by GH here
  # original method assumed that Sigma was already transformed
  Sigma <- diag(stders) %*% cor %*% diag(stders)

  #------------------------
  # Added here from outside

  if( is.null(pkg.env$initialized) ){

    # Load data from file
    # Originally stored in RE2C.RData
    path <- system.file("extdata", package = "remaCor")

    pkg.env$RE2Cor.list <- readRDS(paste(path, 'RE2Cor.list.RDS', sep='/'))
    pkg.env$FEtlow_table <- readRDS(paste(path, 'FEtlow_table.RDS', sep='/'))
    pkg.env$RE2C_table <- readRDS(paste(path, 'RE2C_table.RDS', sep='/'))
    # cal_FEprobs <<- readRDS(paste(path, 'cal_FEprobs.RDS', sep='/'))
    pkg.env$FEprobs <- readRDS(paste(path, 'FEprobs.RDS', sep='/'))
    pkg.env$FEs <- readRDS(paste(path, 'FEs.RDS', sep='/'))
    pkg.env$raw.tau2.zero.prob <- readRDS(paste(path, 'raw.tau2.zero.prob.RDS', sep='/'))
    pkg.env$tau2prob_corlist <- readRDS(paste(path, 'tau2prob_corlist.RDS', sep='/'))
    pkg.env$FEs_ext <- readRDS(paste(path, 'FEs_ext.RDS', sep='/'))  
    pkg.env$initialized <- 1
  }
  pkg.env$initialized = pkg.env$initialized + 1

  ## RE2C uses tabulated matrices to correct stat2_cdf and FEt.lows
  correction.function <- function(correlationMatrix){
    mean_cor <- mean(correlationMatrix[lower.tri(correlationMatrix)])
    int_corr <- floor(mean_cor*10)
    float_corr <- (mean_cor*10) - int_corr
    corr_coeff <- max(1,int_corr+1) # set negitive index values to 1
    aset <- c(corr_coeff,corr_coeff+1)
    #diff(matrix(x)) If x is a matrix then the difference operations are carried out on each column separately.
    RE2C_corr <- RE2Cor[corr_coeff,] + diff(RE2Cor[aset,])*float_corr
    tau2.zero.prob_corr <- tau2prob_cor[corr_coeff]+ diff(tau2prob_cor[aset])*float_corr
    tau2.zero.prob <- pkg.env$raw.tau2.zero.prob*tau2.zero.prob_corr
    return(list(RE2C_corr,tau2.zero.prob))
  }

  nstudy = length(beta)
  
  NR <- test_n( nstudy )

  ## finite sample sizes correction factors
  ##########################################
  RE2Cor <- pkg.env$RE2Cor.list[[min(nstudy-1,6)]]
  tau2prob_cor <- pkg.env$tau2prob_corlist[[min(nstudy-1,6)]]
  
  ## run correction.function(cor)
  correction.list <- correction.function(cor)
  RE2C_corr <- correction.list[[1]];   
  tau2.zero.prob <- correction.list[[2]]

  # set finite sample size correction to nstudy-1 or the max allowed
  nmax = nrow(pkg.env$RE2C_table)
  stat2_cdf = as.numeric(pkg.env$RE2C_table[min(nstudy-1,nmax),])*RE2C_corr

  FEt.lows = as.numeric(pkg.env$FEtlow_table[min(nstudy-1,nmax),])

  isdiag = all(cor[!diag(nrow(cor))] == FALSE)

  # End new
  #-----------------

  # empty return object
   returnValEmpty = c(  stat1         = NA,
                        stat2         = NA,
                        RE2Cp         = NA,
                        RE2Cp.twoStep = NA,
                        QE            = NA,
                        QEp           = NA,
                        Isq           = NA )


  # n study numbers
  # stopifnot (length(beta) == length(stders))
  n = length(beta)
  vars <- stders ** 2
  ws <- 1/vars 
  sigmainv <- tryCatch( solve(Sigma),
    error = function(e){
      warning("At least 1 eigen-value of correlation matrix (cor) is negative,\nso the matrix is not a valid (i.e. positive definite) correlation matrix.\nConsider using Matrix::nearPD().")
      NA
    }
  )  
  if( length(sigmainv) == 1){
    if( is.na(sigmainv) )
      return( t(returnValEmpty ) )
  }
  beta_rv <- matrix(beta,nrow=1,ncol=n)
  ones <- matrix(rep(1,n),nrow=1,ncol=n)
  sumW <- sum(ws)
  sumW2 <- sum(ws ** 2)
  meanBeta <- (ones %*% sigmainv %*% t(beta_rv)) / (ones %*% sigmainv %*% t(ones))
  Q <- (beta_rv - rep(meanBeta,n))** 2 %*% ws
  meanW <- mean(ws)
  Sw2 <- 1/(n-1) * (sumW2 - n*meanW*meanW)
  U = (n-1)*(meanW - Sw2/(n*meanW))
  tau2 <- max( (Q-(n-1))/U, 0 )

  ##-----------------------------------------------
  ## Eigen-decomposition optimization (EMMA- style)
  ##-----------------------------------------------
  K <- Sigma
  eig.L <- eigen(K,symmetric=TRUE)
  L.values <- eig.L$values
  L.vectors <- eig.L$vectors
  S <- diag(n)-matrix(1/n,n,n)
  eig.R <- eigen (S%*%K%*%S,symmetric=TRUE)
  R.values <- eig.R$values[seq(1,n-1)]
  R.vectors <- eig.R$vectors[,seq(1,n-1)]
  etas <- crossprod(R.vectors,t(beta_rv))
  etasqs <- etas^2
  xis <- L.values
  lambdas <- R.values

  if( any(L.values < 0) ){
    warning("At least 1 eigen-value of correlation matrix (cor) is negative,\nso the matrix is not a valid (i.e. positive definite) correlation matrix.\nConsider using Matrix::nearPD().")

      return( t(returnValEmpty ) )
  }else{
    if( any(R.values < 0) ){
      warning("At least 1 eigen-value of correlation matrix (S%*% K %*% S) is negative,\nso the matrix is not a valid (i.e. positive definite) correlation matrix.\nConsider using Matrix::nearPD().")
      return( t(returnValEmpty ) )
    }
  }

  mle.tau2 <- NR(tau2,n,xis,etasqs,lambdas)
  Hinv <- solve(Sigma+mle.tau2*diag(n))
  mle.beta <- (ones %*% Hinv %*% t(beta_rv)) / (ones %*% Hinv %*% t(ones))
  mle.ll <- -LL.fun(mle.tau2,n,xis,etasqs,lambdas)

  null.ll = -ll.fun(c(0,0),beta_rv,Sigma)
  lrt.mle <- -2*(null.ll-mle.ll)

  # get two test statistics
  stat1 = -2*(null.ll+ll.fun(c(meanBeta,0),beta_rv,Sigma))
  stat2 = -2*(-ll.fun(c(meanBeta,0),beta_rv,Sigma) + ll.fun(c(mle.beta,mle.tau2),beta_rv,Sigma))

  if( n <= length(tau2.zero.prob) ){
    # Finite sample size adjustment
    p.RE2_cond <- tau2.zero.prob[n-1]*pchisq(lrt.mle,1,lower.tail=FALSE) + (1-tau2.zero.prob[n-1])*pchisq(lrt.mle,2,lower.tail=FALSE)
  }else{    
    p.RE2_cond <- 0.5*pchisq(lrt.mle,1,lower.tail=FALSE) + 0.5*pchisq(lrt.mle,2,lower.tail=FALSE)
  }
  
  if( ! isdiag ){
    Q = NA
  }

  if( twoStep ){
    ##----------------------------------------------
    ## Given two stats, let's calculate p-value
    ##----------------------------------------------
    stat = stat1 + stat2
    approx_stat1 = min(stat1,49.974)

    RE2C_ext <- function(stat,stat2_cdf,tmax,n){
      extra<-stat-50
      modFEs<-seq(extra+0.25,stat,0.25)
      modFEprobs<-cal_FEprobs(modFEs)
      p.RE2C=(((1-tau2.zero.prob[n-1])*stat2_cdf[floor(20*(apply(cbind(rep(tmax,200),(rep(50,200)-pkg.env$FEs_ext)),1,max)))+1])%*%modFEprobs + (1-tau2.zero.prob[n-1])*stat2_cdf[floor(20*tmax+1)]*pchisq(stat, 1, lower.tail=FALSE))
      return(p.RE2C)
    } 

    if( c(FEt.lows)[floor(approx_stat1*20)+1] > stat2 ){
      p.RE2C.twoStep = 1
      # GEH need to compute this value anyway, right?

    }else if(stat>50){
      p.RE2C.twoStep <- RE2C_ext(stat,stat2_cdf=stat2_cdf,tmax=FEt.lows[1000],n)
        
    }else{   
      p.RE2C.twoStep <- ((1-tau2.zero.prob[n-1])*stat2_cdf[floor(20*(apply(cbind(FEt.lows,(rep(stat,1000)-pkg.env$FEs)),1,max)))+1])%*%pkg.env$FEprobs + (1-tau2.zero.prob[n-1])*stat2_cdf[floor(20*FEt.lows[1000]+1)]*pchisq(50, 1, lower.tail=FALSE)
    }
  }else{
    p.RE2C.twoStep = NA
  }

  returnVal = c(  stat1         = stat1,
                  stat2         = stat2,
                  RE2Cp         = p.RE2_cond,
                  RE2Cp.twoStep = p.RE2C.twoStep,
                  QE            = Q,
                  QEp           = pchisq(q=Q,df=(nstudy-1),lower.tail=FALSE),
                  Isq           = max(100*(Q-(nstudy-1))/Q,0) )

   data.frame( t(returnVal))
}

