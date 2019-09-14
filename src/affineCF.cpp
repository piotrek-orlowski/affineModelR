#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;


arma::cube affineCFevalCpp(const arma::cube coeffs, const arma::mat stateMat, const bool retLog){
  
  // Pull out coeffs dimensions for convenience
  unsigned int Np1, T, U;
  U = coeffs.n_rows;
  T = coeffs.n_cols;
  Np1 = coeffs.n_slices;
  
  // S states
  unsigned int S = stateMat.n_rows;
  
  // extend stateMat
  arma::mat stateMatExtended(S,Np1,arma::fill::ones);
  stateMatExtended.submat(0,0,S-1,Np1-2) = stateMat;
  
  // Every 1x1xNp1 tube = coefficients to obtain log-MGF at maturity t and 
  // argument u. Every such tube has to be multiplied by all rows of stateMat
  // and added together.
  
  // Initialise CF container
  arma::cube cfVals(U,T,S,arma::fill::zeros);
  
  // Initialise local coefficient more container
  arma::mat locCoeffs(U,T,arma::fill::zeros);
  // Loop over states and fill cfVals;
  for(unsigned int tind = 0; tind < T; tind++){
    locCoeffs = coeffs.subcube(0,tind,0,U-1,tind,Np1-1);
    cfVals.subcube(0,tind,0,U-1,tind,S-1) =  locCoeffs * stateMatExtended.t();
  }
  
  if(!retLog){
    cfVals = arma::exp(cfVals); 
  }
  
  return cfVals;
}

// [[Rcpp::export]]
arma::cx_cube affineCFevalCpp(const arma::cx_cube coeffs, const arma::mat stateMat, const bool retLog){
  
  // Pull out coeffs dimensions for convenience
  unsigned int Np1, T, U;
  U = coeffs.n_rows;
  T = coeffs.n_cols;
  Np1 = coeffs.n_slices;
  
  // S states
  unsigned int S = stateMat.n_rows;
  
  // extend stateMat
  arma::mat stateMatExtended(S,Np1,arma::fill::ones);
  stateMatExtended.submat(0,0,S-1,Np1-2) = stateMat;
  
  // Every 1x1xNp1 tube = coefficients to obtain log-MGF at maturity t and 
  // argument u. Every such tube has to be multiplied by all rows of stateMat
  // and added together.
  
  // Initialise CF container
  arma::cx_cube cfVals(U,T,S,arma::fill::zeros);
  
  // Initialise local coefficient more container
  arma::cx_mat locCoeffs(U,T,arma::fill::zeros);
  // Loop over states and fill cfVals;
  for(unsigned int tind = 0; tind < T; tind++){
    locCoeffs = coeffs.subcube(0,tind,0,U-1,tind,Np1-1);
    cfVals.subcube(0,tind,0,U-1,tind,S-1) =  locCoeffs * stateMatExtended.t();
  }
  
  if(!retLog){
    cfVals = arma::exp(cfVals); 
  }
  
  return cfVals;
}