// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "jumpGenerator.h"
#include "../inst/include/affineModelR.h"

//'@export
//[[Rcpp::export]]
SEXP getPointerToGenerator(std::string fstr) {
  if (fstr == "expNormJumpTransform")
    return(Rcpp::wrap(Rcpp::XPtr<funcPtr>(new funcPtr(&generate_expNormJump))));
  else if (fstr == "kouExpJumpTransform")
    return(Rcpp::wrap(Rcpp::XPtr<funcPtr>(new funcPtr(&generate_kouExpJump))));
  else
    return Rcpp::wrap(Rcpp::XPtr<funcPtr>(R_NilValue)); // runtime error as NULL no XPtr
}

//'@export
//[[Rcpp::export]]
arma::vec testPointerToGenerator(SEXP genPtr_, arma::vec jmpPar){
  // Get the jump generator pointer
  Rcpp::XPtr<funcPtr> genPtr(genPtr_);
  funcPtr genFoo = *genPtr;
  arma::vec rndVec = genFoo(jmpPar);
  return(rndVec);
}