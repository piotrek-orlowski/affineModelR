// [[Rcpp::depends(BH)]]

#include "RcppArmadillo.h"
#include "../inst/include/expm1c.h"
#include <boost/math/special_functions/sign.hpp>

using namespace std;
using namespace Rcpp;
//' @export
// [[Rcpp::export]]
arma::vec generate_expNormJump(Rcpp::List jmpPar){
  
  arma::vec res(2,arma::fill::zeros);
  
  arma::vec stockJumpMean(1);
  stockJumpMean(0) = Rcpp::as<double>(jmpPar["muYc"]);
  arma::vec stockJumpVol(1);
  stockJumpVol(0) = Rcpp::as<double>(jmpPar["sigmaYc"]);
  arma::vec volJumpMean(1);
  volJumpMean(0) = 1.0/Rcpp::as<double>(jmpPar["muSc"]);
  arma::vec volJumpCorr(1);
  volJumpCorr(0) = Rcpp::as<double>(jmpPar["rhoc"]);
  
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
arma::vec generate_kouExpJump(Rcpp::List jmpPar){
  
  arma::vec res(2,arma::fill::zeros);
  arma::vec stockJumpMean(1);
  stockJumpMean(0) = Rcpp::as<double>(jmpPar["muYc"]);
  arma::vec stockJumpVol(1);
  stockJumpVol(0) = Rcpp::as<double>(jmpPar["sigmaYc"]);
  arma::vec volJumpMean(1);
  volJumpMean(0) = Rcpp::as<double>(jmpPar["muSc"]);
  arma::vec volJumpCorr(1);
  volJumpCorr(0) = Rcpp::as<double>(jmpPar["rhoc"]);
  
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
arma::vec generate_expStockRetJump(Rcpp::List jmpPar){
  
  // jumps in underlying only (just in case, set a jump of zero in first vol)
  arma::vec res(2,arma::fill::zeros);
  
  // jump parameters
  double stockPar = Rcpp::as<double>(jmpPar["muYc"]);
  int stockParSign = boost::math::sign(stockPar);
  double stockJump = abs(1.0/stockPar);
  
  arma::vec tmpRnd = rexp(1, stockJump) * stockParSign; 
  
  res(0) = tmpRnd(0);
  
  return(res);
}

//' @export
// [[Rcpp::export]]
arma::vec generate_bgl2019Jump(Rcpp::List jmpPar){
  
  // jumps in underlying only (just in case, set a jump of zero in first vol)
  arma::vec res(2,arma::fill::zeros);
  
  // jump parameters
  double muY = Rcpp::as<double>(jmpPar["muY"]);
  double muYinv = 1.0/muY;
  double muV = Rcpp::as<double>(jmpPar["muV"]);
  double muVinv = 1.0/muV;
  
  bool yPresent = (muY != 0.0);
  bool vPresent = (muV != 0.0);
  
  int stockParSign = boost::math::sign(muY);
  int volParSign = boost::math::sign(muV);
  
  muY = stockParSign * muY;
  muV = volParSign * muV;
  
  arma::vec tmpRnd(1);
  tmpRnd.zeros();
  
  if(yPresent){
    tmpRnd = rexp(1, 1.0 / muY) * stockParSign; 
  }
  
  res(0) = tmpRnd(0);
  
  tmpRnd.zeros();
  
  if(vPresent){
    tmpRnd = rexp(1, 1.0 / muV) * volParSign; 
  }
  
  res(1) = tmpRnd(0);
  
  return(res);
}

//' @export
// [[Rcpp::export]]
arma::vec generate_JT2010_cojump(Rcpp::List jmpPar){
  // This is Jacod and Todorov 2010, Do Price and Volatility Jump Together -- co-jumps only
  arma::vec res(2,arma::fill::zeros);
  
  // in jmpPar the first two are l and h, the last two are d and u
  res(0) = Rcpp::as<double>(jmpPar["muYc"]) + Rcpp::as<double>(jmpPar["sigmaYc"]) * arma::randu();
  arma::vec svec(1);
  svec(0) = arma::randu() - 0.5;
  svec = arma::sign(svec);
  res(0) *= svec(0);
  
  res(1) = Rcpp::as<double>(jmpPar["muSc"]) + Rcpp::as<double>(jmpPar["rhoc"]) * arma::randu();
  
  return(res);
}

//' @export
// [[Rcpp::export]]
arma::vec generate_JT2010_cojump_voljump(Rcpp::List jmpPar){
  
  // This is Jacod and Todorov 2010, Do Price and Volatility Jump Together -- co-jumps and pure vol jumps
  arma::vec res(2,arma::fill::zeros);
  
  // check if co-jump or vol jumps
  arma::vec covec(1);
  covec(0) = arma::randu() - 0.5;
  
  // in jmpPar the first two are l and h, the last two are d and u
  if(covec(0) <= 0){
    res(0) = Rcpp::as<double>(jmpPar["muYc"]) + Rcpp::as<double>(jmpPar["sigmaYc"]) * arma::randu();
    arma::vec svec(1);
    svec(0) = arma::randu() - 0.5;
    svec = arma::sign(svec);
    res(0) *= svec(0);
    
    res(1) = Rcpp::as<double>(jmpPar["muSc"]) + Rcpp::as<double>(jmpPar["rhoc"]) * arma::randu();
  } else {
    res(0) = 0.0;
    res(1) = Rcpp::as<double>(jmpPar["muSc"]) + Rcpp::as<double>(jmpPar["rhoc"]) * arma::randu();
  }
  return(res);
}

//' @export
// [[Rcpp::export]]
arma::vec generate_1sidedExp(Rcpp::List jmpPar){
  
  // jumps in underlying and all 3 factors
  arma::vec res(4,arma::fill::zeros);
  
  // jump parameters
  double stockJump = 1.0/Rcpp::as<double>(jmpPar["muStock"]);
  double volJump = 1.0/Rcpp::as<double>(jmpPar["muVol"]);
  double jmpRho = Rcpp::as<double>(jmpPar["rhoc"]);
  double intJump = 1.0/Rcpp::as<double>(jmpPar["muInt"]);
  double vol2Jump = 1.0/Rcpp::as<double>(jmpPar["muVol2"]);
  
  // either asset/1st jump or 2nd/3rd jump
  double gammaProp = Rcpp::as<double>(jmpPar["gammaProp"]);
  bool firstJump = arma::randu() <= 1.0/(1.0 + gammaProp);
  
  if(firstJump){
    double vJmp = rexp(1, volJump)[0];
    double sJmp = jmpRho * vJmp - rexp(1,stockJump)[0];
    res(0) = sJmp;
    res(1) = vJmp;
  } else {
    double iJmp = rexp(1, intJump)[0];
    double v2Jmp = rexp(1,vol2Jump)[0];
    res(2) = iJmp;
    res(3) = v2Jmp;
  }
  
  return(res);
}

