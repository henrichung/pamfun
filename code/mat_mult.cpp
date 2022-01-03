// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP armaMatMult(arma::mat A, arma::mat B){
    arma::mat C = A.t() * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP armaMatMult2(arma::mat A, arma::mat B){
    arma::mat C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP sparseMatMult(const arma::sp_mat A, const arma::sp_mat B){
	arma::sp_mat C = A * B;
	return Rcpp::wrap(C);
}