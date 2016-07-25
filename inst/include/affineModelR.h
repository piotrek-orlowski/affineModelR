#ifndef AFFINE_MODELR_H
#define AFFINE_MODELR_H
#include<RcppArmadillo.h>
// #include "affineModelR_RcppExports.h"

typedef arma::vec (*funcPtr)(const Rcpp::List x);

typedef std::complex<double> (*cmpFuncPtr) (const arma::cx_colvec beta, const Rcpp::List jmpPar);
typedef arma::cx_mat (*cmpFuncPtrMat) (const arma::cx_colvec beta, const Rcpp::List jmpPar);

#endif