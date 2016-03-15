#include <RcppArmadillo.h>
#include "../inst/include/affineModelR.h"
using namespace std;
using namespace Rcpp;

// [[Rcpp::interfaces(r,cpp)]]
// [[Rcpp::export]]
arma::cube affineCFderivsEvalCpp(const arma::cube& coeffs, const arma::mat& stateMat){
  
  // The coeffs cube is of size U x T x 4*(N+1)
  int U = coeffs.n_rows;
  int T = coeffs.n_cols;
  int N = stateMat.n_cols;
  int Np1 = N+1;
  int Nbig = coeffs.n_slices;
  int S = stateMat.n_rows;
  
  // Extract CF coeffs
  arma::cube coeffs_cf = coeffs.subcube(0,0,0,U-1,T-1,Np1-1);
  
  // Calculate cf
  arma::cube cfVals = affineModelR::affineCFevalCpp(coeffs_cf, stateMat, false);
  
  // Extract first derivative coeffs
  arma::cube coeffs_d1 = coeffs.subcube(0,0,1*Np1,U-1,T-1,2*Np1-1);
  
  // calculate derivative
  arma::cube cfD1Vals_log = affineModelR::affineCFevalCpp(coeffs_d1, stateMat, true);
  arma::cube cfD1Vals = cfD1Vals_log % cfVals;
  
  // Extract second derivative coeffs
  arma::cube coeffs_d2 = coeffs.subcube(0,0,2*Np1,U-1,T-1,3*Np1-1);
  
  // calculate second derivative
  arma::cube cfD2Vals_log = affineModelR::affineCFevalCpp(coeffs_d2, stateMat, true);
  arma::cube cfD2Vals = (arma::pow(cfD1Vals_log,2.0) + cfD2Vals_log) % cfVals;
  
  // Extract third derivative coeffs
  arma::cube coeffs_d3 = coeffs.subcube(0,0,3*Np1,U-1,T-1,4*Np1-1);
  
  // calculate third derivative
  arma::cube cfD3Vals_log = affineModelR::affineCFevalCpp(coeffs_d3, stateMat, true);
  arma::cube cfD3Vals = cfVals % (cfD3Vals_log + 3*cfD1Vals_log % cfD2Vals_log + arma::pow(cfD1Vals_log,3.0));
  
  // splice the arrays together
  arma::cube cfRet(U,T,4*S);
  cfRet.slices(0,S-1) = cfVals;
  
  cfRet.slices(S,2*S-1) = cfD1Vals;
  
  cfRet.slices(2*S,3*S-1) = cfD2Vals;
  
  cfRet.slices(3*S,4*S-1) = cfD3Vals;
   
  return cfRet;
}