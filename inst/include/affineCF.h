#ifndef AFFINE_CF_EVAL_H
#define AFFINE_CF_EVAL_H

#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;

arma::cube affineCFevalCpp(const arma::cube coeffs, const arma::mat stateMat, const bool retLog);

arma::cube affineCFderivsEvalCpp(const arma::cube coeffs, const arma::mat stateMat);

arma::cx_cube affineCFevalCpp(const arma::cx_cube coeffs, const arma::mat stateMat, const bool retLog);  

#endif