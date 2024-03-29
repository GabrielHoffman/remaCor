// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rwishart_chol
arma::mat rwishart_chol(const int df, const arma::mat& S_chol);
RcppExport SEXP _remaCor_rwishart_chol(SEXP dfSEXP, SEXP S_cholSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S_chol(S_cholSEXP);
    rcpp_result_gen = Rcpp::wrap(rwishart_chol(df, S_chol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_remaCor_rwishart_chol", (DL_FUNC) &_remaCor_rwishart_chol, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_remaCor(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
