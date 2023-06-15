
#include <RcppArmadillo.h>

// [[Rcpp::export(rng=true)]]
arma::mat rwishart_chol(const int df, const arma::mat& S_chol) {
    arma::uword m = S_chol.n_cols;
    arma::uword i, j;
    arma::mat A(m, m, arma::fill::zeros);
    for ( i = 1; i < m; ++i ) {
        for ( j = 0; j < i; ++j ) {
            A(i, j) = R::rnorm(0.0, 1.0);
        }
    }
    for ( i = 0; i < m; ++i ) {
        A(i, i) = sqrt(R::rchisq(df - i));
    }
    return A.t() * S_chol;
}