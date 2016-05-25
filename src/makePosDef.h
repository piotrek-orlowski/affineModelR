#ifndef _affineCpp_MAKEPOSDEF
#define _affineCpp_MAKEPOSDEF

#include "RcppArmadillo.h"

//This function takes a symmetric and potentially non positive definite matrix and makes it positive definite
void makePosDef(arma::mat& semiDefMat, double relEig = 1e-6);

#endif