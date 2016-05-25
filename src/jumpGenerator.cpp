#include "RcppArmadillo.h"
#include "../inst/include/expm1c.h"

using namespace std;
using namespace Rcpp;
//' @export
// [[Rcpp::export]]
arma::vec generate_expNormJump(arma::vec jmpPar){
  
  arma::vec res(2,arma::fill::zeros);
  
  arma::vec stockJumpMean(1);
  stockJumpMean(0) = jmpPar(0);
  arma::vec stockJumpVol(1);
  stockJumpVol(0) = jmpPar(1);
  arma::vec volJumpMean(1);
  volJumpMean(0) = 1.0/jmpPar(2);
  arma::vec volJumpCorr(1);
  volJumpCorr(0) = jmpPar(3);
  
  arma::vec tmpRnd(1);
  // Generate vol jump
  tmpRnd = rexp(1,volJumpMean(0));
  // tmpRnd = arma::randg(1, arma::distr_param(1.0,volJumpMean(0)));
  res(1) = tmpRnd(0);
  
  // Generate stock jump
  tmpRnd = arma::randn(1) % stockJumpVol;
  tmpRnd(0) += stockJumpMean(0);
  tmpRnd(0) += res(1) * volJumpCorr(0);
  
  res(0) = tmpRnd(0);
  return(res);
}

//' @export
// [[Rcpp::export]]
arma::vec generate_kouExpJump(arma::vec jmpPar){
  
  arma::vec res(2,arma::fill::zeros);
  arma::vec stockJumpMean(1);
  stockJumpMean(0) = jmpPar(0);
  arma::vec stockJumpVol(1);
  stockJumpVol(0) = jmpPar(1);
  arma::vec volJumpMean(1);
  volJumpMean(0) = jmpPar(2);
  arma::vec volJumpCorr(1);
  volJumpCorr(0) = jmpPar(3);
  
  arma::vec tmpRnd(1);
  arma::ivec coinFlip(1);
  // Generate vol jump
  tmpRnd = rexp(1,volJumpMean(0));
  // tmpRnd = arma::randg(1, arma::distr_param(1.0,1.0/volJumpMean(0)));
  // res(1) = as<double>(rgamma(1,1.0,1.0/volJumpMean(0)));
  res(1) = tmpRnd(0);
  
  // Generate coin flip
  coinFlip = arma::randi(1, arma::distr_param(0,1));
  coinFlip(0) = 2 * coinFlip(0) - 1;
  
  // Generate stock jump
  tmpRnd = rexp(1,1.0/stockJumpVol(0));
  // tmpRnd = arma::randg(1, arma::distr_param(1.0,stockJumpVol(0)));
  tmpRnd(0) *= coinFlip(0);
  tmpRnd(0) += stockJumpMean(0);
  tmpRnd(0) += res(1) * volJumpCorr(0);
   
  res(0) = tmpRnd(0);
  
  return(res);
}

//' @export
// [[Rcpp::export]]
arma::vec generate_JT2010_cojump(arma::vec jmpPar){
  // This is Jacod and Todorov 2010, Do Price and Volatility Jump Together -- co-jumps only
  arma::vec res(2,arma::fill::zeros);
  
  // in jmpPar the first two are l and h, the last two are d and u
  res(0) = jmpPar(0) + jmpPar(1) * arma::randu();
  arma::vec svec(1);
  svec(0) = arma::randu() - 0.5;
  svec = arma::sign(svec);
  res(0) *= svec(0);
  
  res(1) = jmpPar(2) + jmpPar(3) * arma::randu();
  
  return(res);
}

//' @export
// [[Rcpp::export]]
arma::vec generate_JT2010_cojump_voljump(arma::vec jmpPar){
  
  // This is Jacod and Todorov 2010, Do Price and Volatility Jump Together -- co-jumps and pure vol jumps
  arma::vec res(2,arma::fill::zeros);
  
  // check if co-jump or vol jumps
  arma::vec covec(1);
  covec(0) = arma::randu() - 0.5;
  
  // in jmpPar the first two are l and h, the last two are d and u
  if(covec(0) <= 0){
    res(0) = jmpPar(0) + jmpPar(1) * arma::randu();
    arma::vec svec(1);
    svec(0) = arma::randu() - 0.5;
    svec = arma::sign(svec);
    res(0) *= svec(0);
    
    res(1) = jmpPar(2) + jmpPar(3) * arma::randu();
  } else {
    res(0) = 0.0;
    res(1) = jmpPar(2) + jmpPar(3) * arma::randu();
  }
  return(res);
}