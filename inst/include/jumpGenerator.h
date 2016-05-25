#ifndef _affineCpp_JGEN
#define _affineCpp_JGEN

#include "RcppArmadillo.h"

arma::vec generate_expNormJump(arma::vec jmpPar);

arma::vec generate_kouExpJump(arma::vec jmpPar);

arma::vec generate_JT2010_cojump(arma::vec jmpPar);

arma::vec generate_JT2010_cojump_voljump(arma::vec jmpPar);

#endif