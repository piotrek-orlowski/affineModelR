#ifndef _affineCpp_JTR
#define _affineCpp_JTR

#include "RcppArmadillo.h"
#include <complex>

std::complex<double> jumpTransform(const arma::cx_colvec& beta, const double& muYc, const double& sigmaYc, const double& muSc, const double& rhoc);

arma::cx_mat jumpTransformD1(const arma::cx_colvec& beta, const double& muYc, const double& sigmaYc, const double& muSc, const double& rhoc);

arma::cx_mat jumpTransformD2(const arma::cx_colvec& beta, const double& muYc, const double& sigmaYc, const double& muSc, const double& rhoc);

arma::cx_mat jumpTransformD3(const arma::cx_colvec& beta, const double& muYc, const double& sigmaYc, const double& muSc, const double& rhoc);

std::complex<double> kouExpTransform(const arma::cx_colvec& beta, const double& muYc, const double& sigmaYc, const double& muSc, const double& rhoc);

arma::cx_mat kouExpTransformD1(const arma::cx_colvec& beta, const double& muYc, const double& sigmaYc, const double& muSc, const double& rhoc);

arma::cx_mat kouExpTransformD2(const arma::cx_colvec& beta, const double& muYc, const double& sigmaYc, const double& muSc, const double& rhoc);

arma::cx_mat kouExpTransformD3(const arma::cx_colvec& beta, const double& muYc, const double& sigmaYc, const double& muSc, const double& rhoc);
#endif