#ifndef _affineCpp_JTR
#define _affineCpp_JTR

#include "RcppArmadillo.h"
#include <complex>

std::complex<double> jumpTransform(const arma::cx_colvec& beta, const Rcpp::List& jmpPar);

arma::cx_mat jumpTransformD1(const arma::cx_colvec& beta, const Rcpp::List& jmpPar);

arma::cx_mat jumpTransformD2(const arma::cx_colvec& beta, const Rcpp::List& jmpPar);

arma::cx_mat jumpTransformD3(const arma::cx_colvec& beta, const Rcpp::List& jmpPar);

std::complex<double> kouExpTransform(const arma::cx_colvec& beta, const Rcpp::List& jmpPar);

arma::cx_mat kouExpTransformD1(const arma::cx_colvec& beta, const Rcpp::List& jmpPar);

arma::cx_mat kouExpTransformD2(const arma::cx_colvec& beta, const Rcpp::List& jmpPar);

arma::cx_mat kouExpTransformD3(const arma::cx_colvec& beta, const Rcpp::List& jmpPar);
#endif