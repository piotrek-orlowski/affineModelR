#ifndef PASS_H
#define PASS_H

#include<RcppArmadillo.h>

typedef arma::vec (*funcPtr)(const arma::vec x);

typedef std::complex<double> (*cmpFuncPtr) (const arma::cx_colvec& beta, const double& muYc, const double& sigmaYc, const double& muSc, const double& rhoc);
typedef arma::cx_mat (*cmpFuncPtrMat) (const arma::cx_colvec& beta, const double& muYc, const double& sigmaYc, const double& muSc, const double& rhoc);

#endif