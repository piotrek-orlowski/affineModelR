#ifndef AFFINE_CF_EVAL_H
#define AFFINE_CF_EVAL_H

#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;

arma::cube affineCFevalCpp(const NumericVector coeffs, const arma::mat stateMat);

arma::cube affineCFderivsEvalCpp(const arma::cube coeffs, const arma::mat stateMat);
  
// arma::cx_cube affineCFevalCpp(const arma::cx_cube& coeffs, const arma::mat& stateMat);
// 
// arma::cube affineCFderivsEvalCpp(const arma::cube& coeffs, const arma::mat& stateMat);
// 
// arma::cx_cube affineCFDerivsEvalCpp(const arma::cx_cube& coeffs, const arma::mat& stateMat);

#endif