
library(RUnit)
library(metafor)

test_fixed_effect = function(){

	set.seed(1)
	n = 1000
	beta = rnorm(n)
	se = runif(n, .1, 1)

	# Lin Sullivan meta-analysis with identity correlation
	resLS = LS( beta, se)

	# standard fixed effects meta analysis
	resRMA = rma( yi=beta, sei=se, method="FE" )

	# check that results are the same
	checkEquals(resLS$beta, as.numeric(resRMA$beta) ) & checkEquals(resLS$se, as.numeric(resRMA$se) ) & checkEquals(resLS$p, resRMA$pval)
}


test_random_effect = function(){

	set.seed(1)
	n = 49
	beta = rnorm(n)
	se = runif(n, .1, 1)

	# Lin Sullivan meta-analysis with identity correlation
	resRE2C = RE2C( beta, se)

	# standard random effects meta analysis
	resRMA = rma( yi=beta, sei=se, method="REML" )

	# check that results are the same
	checkEquals(resRE2C$QE, resRMA$QE ) & checkEquals(resRE2C$QEp, resRMA$QEp ) 
}

