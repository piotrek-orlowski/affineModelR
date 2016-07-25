#ifndef _affineCpp_JGEN
#define _affineCpp_JGEN

#include "RcppArmadillo.h"

arma::vec generate_expNormJump(Rcpp::List jmpPar);

arma::vec generate_kouExpJump(Rcpp::List jmpPar);

arma::vec generate_JT2010_cojump(Rcpp::List jmpPar);

arma::vec generate_JT2010_cojump_voljump(Rcpp::List jmpPar);

arma::vec generate_1sidedExp(Rcpp::List jmpPar);

#endif