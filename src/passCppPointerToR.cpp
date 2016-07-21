// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "../inst/include/jumpGenerator.h"
#include "../inst/include/jumpTransform.h"
#include "../inst/include/affineModelR.h"

//'@export
//[[Rcpp::export]]
SEXP getPointerToGenerator(std::string fstr) {
  if (fstr == "expNormJumpTransform")
    return(Rcpp::wrap(Rcpp::XPtr<funcPtr>(new funcPtr(&generate_expNormJump))));
  else if (fstr == "kouExpJumpTransform")
    return(Rcpp::wrap(Rcpp::XPtr<funcPtr>(new funcPtr(&generate_kouExpJump))));
  else if (fstr == "jacodTodorovCoJumps")
    return(Rcpp::wrap(Rcpp::XPtr<funcPtr>(new funcPtr(&generate_JT2010_cojump))));
  else if (fstr == "jacodTodorovCoJumpsVolJumps")
    return(Rcpp::wrap(Rcpp::XPtr<funcPtr>(new funcPtr(&generate_JT2010_cojump_voljump))));
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
  else if (fstr == "jacodTodorovCoJumps"){
    pointers = Rcpp::List::create(
      Rcpp::Named("TF") = Rcpp::XPtr<cmpFuncPtr>(new cmpFuncPtr(&jt2010_transform_CJ))
    );
    return(pointers);
  }
  else if (fstr == "jacodTodorovCoJumpsVolJumps"){
    pointers = Rcpp::List::create(
      Rcpp::Named("TF") = Rcpp::XPtr<cmpFuncPtr>(new cmpFuncPtr(&jt2010_transform_CJ_VJ))
    );
    return(pointers);
  }
  else
    return Rcpp::wrap(Rcpp::XPtr<cmpFuncPtr>(R_NilValue)); // runtime error as NULL no XPtr
}


//'@export
//[[Rcpp::export]]
arma::vec evaluateGenerator(SEXP genPtr_, arma::vec jmpPar){
  // Get the jump generator pointer
  Rcpp::XPtr<funcPtr> genPtr(genPtr_);
  funcPtr genFoo = *genPtr;
  arma::vec rndVec = genFoo(jmpPar);
  return(rndVec);
}

//'@export
//[[Rcpp::export]]
arma::cx_colvec evaluateTransform(SEXP genPtr_, const arma::cx_colvec beta, const Rcpp::List jmpPar){
// std::complex<double> evaluateTransform(SEXP genPtr_, const arma::cx_colvec beta, const Rcpp::List jmpPar){
  // Get the jump generator pointer
  Rcpp::XPtr<cmpFuncPtr> genPtr(genPtr_);
  cmpFuncPtr genFoo = *genPtr;
  // std::complex<double> tfVal = genFoo(beta,jmpPar);
  arma::cx_colvec tfVal = genFoo(beta,jmpPar);
  return(tfVal);
}


// //'@export
// //[[Rcpp::export]]
// arma::cx_mat evaluateTransformDerivative(SEXP genPtr_, const arma::cx_colvec& beta, const Rcpp::List& jmpPar){
//   // Get the jump generator pointer
//   Rcpp::XPtr<cmpFuncPtrMat> genPtr(genPtr_);
//   cmpFuncPtr genFoo = *genPtr;
//   arma::cx_mat tfVal = genFoo(beta,jmpPar);
//   return(tfVal);
// }