// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace RcppEigen;

// [[Rcpp::export]]
SEXP e(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B, Eigen::Map<Eigen::MatrixXd> C){
  Eigen::setNbThreads(1);
  Eigen::MatrixXd D = A * B * C;

  return Rcpp::wrap(D);
}

