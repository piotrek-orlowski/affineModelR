// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "jumpGenerator.h"
#include "jumpTransform.h"
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
SEXP getPointerToJumpTransform(std::string fstr) {
  Rcpp::List pointers;
  if (fstr == "expNormJumpTransform"){
    pointers = Rcpp::List::create(
      Rcpp::Named("TF") = Rcpp::XPtr<cmpFuncPtr>(new cmpFuncPtr(&jumpTransform)), 
      Rcpp::Named("D1") = Rcpp::XPtr<cmpFuncPtrMat>(new cmpFuncPtrMat(&jumpTransformD1)),
      Rcpp::Named("D2") = Rcpp::XPtr<cmpFuncPtrMat>(new cmpFuncPtrMat(&jumpTransformD2)),
      Rcpp::Named("D3") = Rcpp::XPtr<cmpFuncPtrMat>(new cmpFuncPtrMat(&jumpTransformD3))
    );
    return(pointers);
  }
  else if (fstr == "kouExpJumpTransform"){
    pointers = Rcpp::List::create(
      Rcpp::Named("TF") = Rcpp::XPtr<cmpFuncPtr>(new cmpFuncPtr(&kouExpTransform)), 
      Rcpp::Named("D1") = Rcpp::XPtr<cmpFuncPtrMat>(new cmpFuncPtrMat(&kouExpTransformD1)),
      Rcpp::Named("D2") = Rcpp::XPtr<cmpFuncPtrMat>(new cmpFuncPtrMat(&kouExpTransformD2)),
      Rcpp::Named("D3") = Rcpp::XPtr<cmpFuncPtrMat>(new cmpFuncPtrMat(&kouExpTransformD3))
    );
    return(pointers);
  }
  else
    return Rcpp::wrap(Rcpp::XPtr<cmpFuncPtr>(R_NilValue)); // runtime error as NULL no XPtr
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