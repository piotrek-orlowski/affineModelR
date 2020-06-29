// Taylor expansion of the jump transform
#include "RcppArmadillo.h"
#include "../inst/include/expm1c.h"

using namespace std;
// [[Rcpp::export]]
std::complex<double> jumpTransform(const arma::cx_colvec beta, const Rcpp::List jmpPar) {
  try{
    double muYc = Rcpp::as<double>(jmpPar["muYc"]);
    double sigmaYc = Rcpp::as<double>(jmpPar["sigmaYc"]);
    double rhoc = Rcpp::as<double>(jmpPar["rhoc"]);
    double muSc = Rcpp::as<double>(jmpPar["muSc"]);
    // assumption is that only the underlying and the first volatility factor can jumps
    complex<double> c = muSc * beta(1) + rhoc * muSc * beta(0);
    
    complex<double> cf = expm1c(muYc * beta(0) + 0.5 * pow(sigmaYc * beta(0),2));
    cf += c;
    cf = cf / (complex<double>(1,0) - c);
    
    return(cf);
    
  } catch( std::exception& __ex__) {
    forward_exception_to_r(__ex__);
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return 0;
}

// This returns the gradient of the jumpTransform
// [[Rcpp::export]]
arma::cx_mat jumpTransformD1(const arma::cx_colvec beta, const Rcpp::List jmpPar) {
  try{
    double muYc = Rcpp::as<double>(jmpPar["muYc"]);
    double sigmaYc = Rcpp::as<double>(jmpPar["sigmaYc"]);
    double rhoc = Rcpp::as<double>(jmpPar["rhoc"]);
    double muSc = Rcpp::as<double>(jmpPar["muSc"]);
    
    // assumption is that only the underlying and the first volatility factor can jumps
    complex<double> tmpVar;
    int mDim = beta.size();
    arma::cx_mat jtJacobian(1,mDim,arma::fill::zeros);
    
    // deriv wrt beta(0)
    tmpVar = (1.0 + expm1c(muYc * beta(0) + 0.5 * pow(sigmaYc * beta(0),2)));
    tmpVar *= (muSc * rhoc - muYc*(beta(0)*muSc*rhoc+beta(1)*muSc-1.0) - beta(0)*pow(sigmaYc,2)*(beta(0)*muSc*rhoc+beta(1)*muSc-1.0));
    tmpVar /= pow(beta(0)*muSc*rhoc+beta(1)*muSc-1.0,2);
    
    jtJacobian(0,0) = tmpVar;
    
    tmpVar = 0;
    
    // deriv wrt beta(1)
    tmpVar = muSc * (1.0 + expm1c(muYc * beta(0) + 0.5 * pow(sigmaYc * beta(0),2)));
    tmpVar /= pow(beta(0)*muSc*rhoc+beta(1)*muSc-1.0,2);
    
    jtJacobian(0,1) = tmpVar;
    
    // All subsequent derivatives are 0
    
    return(jtJacobian);
    
  } catch( std::exception& __ex__) {
    forward_exception_to_r(__ex__);
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return 0;
}

// This returns the Hessian of the jumpTransform
// [[Rcpp::export]]
arma::cx_mat jumpTransformD2(const arma::cx_colvec beta, const Rcpp::List jmpPar) {
  try{
    double muYc = Rcpp::as<double>(jmpPar["muYc"]);
    double sigmaYc = Rcpp::as<double>(jmpPar["sigmaYc"]);
    double rhoc = Rcpp::as<double>(jmpPar["rhoc"]);
    double muSc = Rcpp::as<double>(jmpPar["muSc"]);
    // assumption is that only the underlying and the first volatility factor can jumps
    complex<double> tmpVar;
    int mDim = beta.size();
    arma::cx_mat jtHessian(mDim,mDim,arma::fill::zeros);
    
    // 2nd deriv wrt beta(0)
    complex<double> tmp1, tmp2, tmp3;
    tmp1 = 2.0 * pow(muSc * rhoc,2.0) * exp(beta(0)*muYc + 0.5*pow(sigmaYc * beta(0),2.0))/pow(1.0 - beta(1) * muSc - beta(0) * rhoc * muSc,3.0);
    tmp2 = 2.0 * exp(beta(0)*muYc + 0.5*pow(sigmaYc * beta(0),2.0)) * rhoc * muSc * (muYc + beta(0)*pow(sigmaYc,2.0))/pow(1.0 - beta(1) * muSc - beta(0) * rhoc * muSc,2.0);
    tmp3 = (pow(sigmaYc,2.0)*exp(beta(0)*muYc+0.5*pow(sigmaYc*beta(0),2.0))+exp(beta(0)*muYc+0.5*pow(sigmaYc*beta(0),2.0))*pow(muYc + beta(0)*pow(sigmaYc,2.0),2.0))/(1.0 - beta(1) * muSc - beta(0) * rhoc * muSc);
    tmpVar = tmp1 + tmp2 + tmp3;
  
    jtHessian(0,0) = tmpVar;
    
    tmpVar = 0;
    
    // 2nd deriv wrt beta(1)
    tmpVar = (1.0 + expm1c(muYc * beta(0) + 0.5 * pow(sigmaYc * beta(0),2)))*2.0*pow(muSc,2)/pow(-beta(0)*muSc*rhoc -beta(1)*muSc + 1.0,3);
    
    jtHessian(1,1) = tmpVar;
    
    tmpVar = 0;
    
    // mixed derivative
    tmpVar = 2.0*pow(muSc,2)*rhoc*(1.0 + expm1c(muYc * beta(0) + 0.5 * pow(sigmaYc * beta(0),2)))/pow(-beta(0)*muSc*rhoc -beta(1)*muSc + 1.0,3);
    tmpVar += (1.0 + expm1c(muYc * beta(0) + 0.5 * pow(sigmaYc * beta(0),2)))*muSc*(muYc + beta(0)*pow(sigmaYc,2))/pow(-beta(0)*muSc*rhoc -beta(1)*muSc + 1.0,2);
    
    jtHessian(1,0) = tmpVar;
    jtHessian(0,1) = tmpVar;
    // All subsequent derivatives are 0
    
    return(jtHessian);
    
  } catch( std::exception& __ex__) {
    forward_exception_to_r(__ex__);
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return 0;
}


// This returns the derivative of the vectorized Hessian of the jumpTransform, i.e. the third deriv
// [[Rcpp::export]]
arma::cx_mat jumpTransformD3(const arma::cx_colvec beta, const Rcpp::List jmpPar) {
  try{
    double muYc = Rcpp::as<double>(jmpPar["muYc"]);
    double sigmaYc = Rcpp::as<double>(jmpPar["sigmaYc"]);
    double rhoc = Rcpp::as<double>(jmpPar["rhoc"]);
    double muSc = Rcpp::as<double>(jmpPar["muSc"]);
    // assumption is that only the underlying and the first volatility factor can jumps
    complex<double> tmpVar;
    int mDim = beta.size();
    arma::cx_mat jtD3(mDim*mDim,mDim,arma::fill::zeros);
    
    // commonly ocurring variables
    complex<double> expb0 = (1.0 + expm1c(muYc * beta(0) + 0.5 * pow(sigmaYc * beta(0),2)));
    complex<double> denom = 1.0 - beta(0) * muSc * rhoc - beta(1) * muSc;
    complex<double> mvr = muSc * rhoc;
    complex<double> b0sq = muYc + beta(0)*pow(sigmaYc,2);
    
    // 3rd deriv wrt beta(0)
    tmpVar = 6.0*expb0*pow(mvr,3)/pow(denom,4);
    tmpVar += 6.0*expb0*pow(mvr,2)*b0sq/pow(denom,3);
    tmpVar += 3.0*mvr*(expb0*pow(sigmaYc,2) + expb0*pow(b0sq,2))/pow(denom,2);
    tmpVar += expb0*(3.0*pow(sigmaYc,2)*b0sq + pow(b0sq,3))/denom;
    
    jtD3(0,0) = tmpVar;
    
    tmpVar = 0;
    
    // 3rd deriv wrt beta(1)
    tmpVar = 6.0*expb0*pow(muSc,3)/pow(denom,4);
    
    jtD3(mDim+1,1) = tmpVar;
    
    tmpVar = 0;
    
    // 2nd deriv wrt beta(0), then 1st wrt beta(1)
    tmpVar = 6.0*expb0*pow(mvr,2)*muSc/pow(denom,4);
    tmpVar += 4.0*expb0*mvr*muSc*b0sq/pow(denom,3);
    tmpVar += muSc*expb0*(pow(sigmaYc,2)+ pow(b0sq,2))/pow(denom,2);
    
    jtD3(0,1) = tmpVar;
    jtD3(1,0) = tmpVar;
    jtD3(mDim,0) = tmpVar;
    
    tmpVar = 0;
    
    // 2nd deriv wrt beta(1), then 1st wrt beta(0)
    tmpVar = 6.0*expb0*mvr*pow(muSc,2)/pow(denom,4);
    tmpVar += 2.0*expb0*pow(muSc,2)*b0sq/pow(denom,3);
    
    jtD3(1,1) = tmpVar;
    jtD3(mDim,1) = tmpVar;
    jtD3(mDim+1,0) = tmpVar;
    // All subsequent derivatives are 0
    
    return(jtD3);
    
  } catch( std::exception& __ex__) {
    forward_exception_to_r(__ex__);
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return 0;
}

// Kou-Exponential jump transform
// [[Rcpp::export]]
std::complex<double> kouExpTransform(const arma::cx_colvec beta, const Rcpp::List jmpPar){
  try{
    double muYc = Rcpp::as<double>(jmpPar["muYc"]);
    double sigmaYc = Rcpp::as<double>(jmpPar["sigmaYc"]);
    double rhoc = Rcpp::as<double>(jmpPar["rhoc"]);
    double muSc = Rcpp::as<double>(jmpPar["muSc"]);
    // assumption is that only the underlying and the first volatility factor can jump
    complex<double> c = (pow(sigmaYc * beta(0),2.0) - complex<double>(1,0)) * (beta(1) - muSc + rhoc*beta(0));
    
    complex<double> cf = muSc * exp(beta(0) * muYc);
    cf = cf / c - complex<double>(1.0,0.0);
    
    return(cf);
    
  } catch( std::exception& __ex__) {
    forward_exception_to_r(__ex__);
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return 0;
}

// This returns the gradient of kouExpTransform
// [[Rcpp::export]]
arma::cx_mat kouExpTransformD1(const arma::cx_colvec beta, const Rcpp::List jmpPar) {
  try{
    double muYc = Rcpp::as<double>(jmpPar["muYc"]);
    double sigmaYc = Rcpp::as<double>(jmpPar["sigmaYc"]);
    double rhoc = Rcpp::as<double>(jmpPar["rhoc"]);
    double muSc = Rcpp::as<double>(jmpPar["muSc"]);
    // assumption is that only the underlying and the first volatility factor can jumps
    complex<double> tmpVar1, tmpVar2, tmpVar3;
    int mDim = beta.size();
    arma::cx_mat jtJacobian(1,mDim,arma::fill::zeros);
    
    // deriv wrt beta(0)
    tmpVar1 = -muSc * rhoc * exp(beta(0)*muYc);
    tmpVar1 /= (pow(beta(0)*sigmaYc,2.0)-complex<double>(1,0))*pow((beta(1) - muSc + rhoc * beta(0)),2.0);
    tmpVar2 = -2.0 * pow(sigmaYc,2.0) * muSc * beta(0) * exp(beta(0) * muYc);
    tmpVar2 /= pow((pow(sigmaYc*beta(0),2.0) - complex<double>(1,0)),2.0)*(beta(1)-muSc + beta(0)*rhoc);
    tmpVar3 = muSc * muYc * exp(beta(0) * muYc);
    tmpVar3 /= pow((pow(sigmaYc*beta(0),2.0) - complex<double>(1,0)),1.0)*(beta(1)-muSc + beta(0)*rhoc);
    
    jtJacobian(0,0) = tmpVar1 + tmpVar2 + tmpVar3;
    
    tmpVar1 = 0;
    tmpVar2 = 0;
    
    // deriv wrt beta(1)
    tmpVar1 = - muSc * exp(muYc * beta(0));
    tmpVar1 /= (pow(beta(0)*sigmaYc,2.0) - complex<double>(1,0))*pow(beta(1) - muSc + beta(0)*rhoc,2.0);
    
    jtJacobian(0,1) = tmpVar1;
    
    // All subsequent derivatives are 0
    
    return(jtJacobian);
    
  } catch( std::exception& __ex__) {
    forward_exception_to_r(__ex__);
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return 0;
}

// This returns the Hessian of kouExpTransform
// [[Rcpp::export]]
arma::cx_mat kouExpTransformD2(const arma::cx_colvec beta, const Rcpp::List jmpPar) {
  try{
    double muYc = Rcpp::as<double>(jmpPar["muYc"]);
    double sigmaYc = Rcpp::as<double>(jmpPar["sigmaYc"]);
    double rhoc = Rcpp::as<double>(jmpPar["rhoc"]);
    double muSc = Rcpp::as<double>(jmpPar["muSc"]);
    // assumption is that only the underlying and the first volatility factor can jumps
    complex<double> tmpVar1, tmpVar2, tmpVar3;
    int mDim = beta.size();
    arma::cx_mat jtHessian(mDim,mDim,arma::fill::zeros);
    
    // 2nd deriv wrt beta(0)
    tmpVar1 = pow(beta(0)*sigmaYc,2.0) - complex<double>(1,0);
    tmpVar2 = beta(1) - muSc + beta(0)*rhoc;
    tmpVar3 = 8.0 * pow(beta(0),2.0) * pow(sigmaYc,4.0) / pow(tmpVar1,3.0);
    tmpVar3 -= 2*pow(sigmaYc,2.0)/pow(tmpVar1,2.0);
    tmpVar3 *= exp(beta(0) * muYc);
    tmpVar3 -= 4.0 * pow(sigmaYc,2.0) * beta(0) * muYc * exp(muYc * beta(0))/pow(tmpVar1,2.0);
    tmpVar3 += pow(muYc,2.0) * exp(muYc * beta(0))/tmpVar1;
    tmpVar3 *= muSc / tmpVar2;
    tmpVar3 -= 2*muSc*rhoc*(-2.0*pow(sigmaYc,2.0)*beta(0)*exp(muYc*beta(0))/pow(tmpVar1,2.0) + muYc*exp(muYc * beta(0))/tmpVar1)/ (pow(tmpVar2,2.0));
    tmpVar3 += 2.0 * muSc * pow(rhoc,2.0) * exp(beta(0) * muYc) / (tmpVar1 * pow(tmpVar2,3.0));
    
    jtHessian(0,0) = tmpVar3;
    
    tmpVar3 = 0;
    
    // 2nd deriv wrt beta(1)
    tmpVar3 = 2.0 * muSc * exp(muYc * beta(0)) / (tmpVar1 * pow(tmpVar2,3.0));
    
    jtHessian(1,1) = tmpVar3;
    
    tmpVar3 = 0;
    
    // mixed derivative
    tmpVar3 = - muSc * muYc * exp(beta(0) * muYc) / (tmpVar1 * pow(tmpVar2,2.0));
    tmpVar3 += 2.0 * pow(sigmaYc,2.0) * beta(0) * muSc * exp(muYc * beta(0)) / (pow(tmpVar1,2.0)*pow(tmpVar2,2.0));
    tmpVar3 += 2.0 * exp(beta(0)*muYc) * muSc * rhoc / (tmpVar1 * pow(tmpVar2,3.0));
    
    jtHessian(1,0) = tmpVar3;
    jtHessian(0,1) = tmpVar3;
    // All subsequent derivatives are 0
    
    return(jtHessian);
    
  } catch( std::exception& __ex__) {
    forward_exception_to_r(__ex__);
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return 0;
}

// This returns the derivative of the vectorized Hessian of the kouExpTransform, i.e. the third deriv
// [[Rcpp::export]]
arma::cx_mat kouExpTransformD3(const arma::cx_colvec beta, const Rcpp::List jmpPar) {
  try{
    double muYc = Rcpp::as<double>(jmpPar["muYc"]);
    double sigmaYc = Rcpp::as<double>(jmpPar["sigmaYc"]);
    double rhoc = Rcpp::as<double>(jmpPar["rhoc"]);
    double muSc = Rcpp::as<double>(jmpPar["muSc"]);
    // assumption is that only the underlying and the first volatility factor can jump
    complex<double> tmpVar1, tmpVar2, tmpVar3;
    int mDim = beta.size();
    arma::cx_mat jtD3(mDim*mDim,mDim,arma::fill::zeros);
    
    // commonly ocurring variables
    tmpVar1 = pow(beta(0)*sigmaYc,2.0) - complex<double>(1,0);
    tmpVar2 = beta(1) - muSc + beta(0)*rhoc;
    
    // 3rd deriv wrt beta(0)
    tmpVar3 = (-48.0 * pow(sigmaYc,6.0) * pow(beta(0),3.0)/pow(tmpVar1,4.0) + 24.0*pow(sigmaYc,4.0)*beta(0)/pow(tmpVar1,3.0)) * exp(beta(0)*muYc);
    tmpVar3 += 3.0 * muYc * exp(beta(0) * muYc) * (8.0*pow(sigmaYc,4.0)*pow(beta(0),2.0)/pow(tmpVar1,3.0) - 2.0*pow(sigmaYc,2.0)/pow(tmpVar1,2.0));
    tmpVar3 += -6.0 * pow(sigmaYc,2.0) * beta(0) * exp(muYc * beta(0)) * pow(muYc,2.0)/pow(tmpVar1,2.0) + exp(muYc * beta(0)) * pow(muYc,3.0)/tmpVar1;
    tmpVar3 *= muSc/tmpVar2;
    tmpVar3 -= 3.0 * muSc * rhoc *((8.0*pow(sigmaYc,4.0)*pow(beta(0),2.0)/pow(tmpVar1,3.0) - 2 * pow(sigmaYc,2.0)/pow(tmpVar1,2.0))*exp(beta(0)*muYc) - 4.0*pow(sigmaYc,2.0)*beta(0)*exp(muYc*beta(0))*muYc/pow(tmpVar1,2.0) + exp(beta(0)*muYc)*pow(muYc,2.0)/tmpVar1)/pow(tmpVar2,2.0);
    tmpVar3 += 6.0 * muSc * pow(rhoc,2.0) * (-2.0*pow(sigmaYc,2.0)*beta(0)*exp(beta(0)*muYc)/pow(tmpVar1,2.0) + muYc*exp(muYc*beta(0))/tmpVar1)/pow(tmpVar2,3.0);
    tmpVar3 -= 6.0 * muSc * pow(rhoc,3.0) * exp(muYc * beta(0)) / (tmpVar1 * pow(tmpVar2,4.0));
    
    jtD3(0,0) = tmpVar3;
    
    tmpVar3 = 0;
    
    // 3rd deriv wrt beta(1)
    tmpVar3 = -6.0 * muSc * exp(beta(0)*muYc)/(tmpVar1 * pow(tmpVar2,4.0));
    
    jtD3(mDim+1,1) = tmpVar3;
    
    tmpVar3 = 0;
    
    // 2nd deriv wrt beta(0), then 1st wrt beta(1)
    tmpVar3 = - muSc*((8.0 * pow(sigmaYc,4.0) * pow(beta(0),2.0)/pow(tmpVar1,3.0) - 2*pow(sigmaYc,2.0)/pow(tmpVar1,2.0))*exp(beta(0)*muYc) - 4.0*pow(sigmaYc,2.0)*beta(0)*exp(beta(0)*muYc)*muYc/pow(tmpVar1,2.0) + exp(beta(0)*muYc)*pow(muYc,2.0)/tmpVar1)/pow(tmpVar2,2.0);
    tmpVar3 += 4.0 * muSc * rhoc *( -2.0*pow(sigmaYc,2.0)*beta(0)*exp(beta(0)*muYc)/pow(tmpVar1,2.0) + exp(muYc*beta(0))*muYc/tmpVar1) / pow(tmpVar2,3.0);
    tmpVar3 -= 6.0 * muSc * pow(rhoc,2.0) *exp(muYc*beta(0)) / (tmpVar1 * pow(tmpVar2, 4.0));
    
    jtD3(0,1) = tmpVar3;
    jtD3(1,0) = tmpVar3;
    jtD3(mDim,0) = tmpVar3;
    
    tmpVar3 = 0;
    
    // 2nd deriv wrt beta(1), then 1st wrt beta(0)
    tmpVar3 = 2.0 * exp(muYc * beta(0)) * muSc * muYc / (tmpVar1 * pow(tmpVar2,3.0));
    tmpVar3 -= 4.0 * pow(sigmaYc,2.0) * beta(0) * muSc *exp(muYc*beta(0)) / (pow(tmpVar1,2.0)*pow(tmpVar2,3.0));
    tmpVar3 -= 6.0 * muSc * rhoc * exp(beta(0)*muYc) /(tmpVar1 * pow(tmpVar2,4.0));
    
    jtD3(1,1) = tmpVar3;
    jtD3(mDim,1) = tmpVar3;
    jtD3(mDim+1,0) = tmpVar3;
    // All subsequent derivatives are 0
    
    return(jtD3);
    
  } catch( std::exception& __ex__) {
    forward_exception_to_r(__ex__);
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return 0;
}

// Jacod and Todorov (2010): doubly uniform jumps, co-jump only
// [[Rcpp::export]]
std::complex<double> jt2010_transform_CJ(const arma::cx_colvec beta, const Rcpp::List jmpPar){
  
  double l,h,d,u;
  l = Rcpp::as<double>(jmpPar["muYc"]);
  h = Rcpp::as<double>(jmpPar["sigmaYc"]);
  d = Rcpp::as<double>(jmpPar["muSc"]);
  u = Rcpp::as<double>(jmpPar["rhoc"]);
  
  complex<double> c1, c2;
  c1 = beta(0);
  c2 = beta(1);
  
  
  complex<double> cf;
  if(h > 0 & u > 0){
    cf = 1.0/(2.0*(h-l)*(u-d)); 
  } else if(h > 0 & u == 0){
    cf = 1.0/(2.0*(h-l)); 
  } else if(h == 0 & u > 0){
    cf = 1.0/(u-d); 
  } else {
    cf = 1.0;
  }
  
  if(u > 0){
    if(abs(c2) > pow(std::numeric_limits<double>::epsilon(), 0.5)){
      cf *= 1.0/c2*(exp(c2 * u) - exp(c2 * d)); 
    } else {
      cf *= (u-d);
    }  
  }
  
  if(h > 0){
    if(abs(c1) > pow(std::numeric_limits<double>::epsilon(), 0.5)){
      cf *= 1.0/c1*(exp(c1 * h) - exp(c1 * l) + exp(-c1*l) - exp(-c1*h));
    } else {
      cf *= 2*(h-l);
    }  
  }
  
  cf -= 1.0;
  
  return cf;
}

// Jacod and Todorov (2010): doubly uniform jumps, co-jump and vol jump
// [[Rcpp::export]]
std::complex<double> jt2010_transform_CJ_VJ(const arma::cx_colvec beta, const Rcpp::List jmpPar){
  
  double l,h,d,u;
  l = Rcpp::as<double>(jmpPar["muYc"]);
  h = Rcpp::as<double>(jmpPar["sigmaYc"]);
  d = Rcpp::as<double>(jmpPar["muSc"]);
  u = Rcpp::as<double>(jmpPar["rhoc"]);
  
  // First do the co-jump CF
  complex<double> c1, c2;
  c1 = beta(0);
  c2 = beta(1);
  
  complex<double> cf;
  if(h > 0 & u > 0){
    cf = 1.0/(2.0*(h-l)*(u-d)); 
  } else if(h > 0 & u == 0){
    cf = 1.0/(2.0*(h-l)); 
  } else if(h == 0 & u > 0){
    cf = 1.0/(u-d); 
  } else {
    cf = 1.0;
  }
  
  if(u > 0){
    if(abs(c2) > pow(std::numeric_limits<double>::epsilon(), 0.5)){
      cf *= 1.0/c2*(exp(c2 * u) - exp(c2 * d)); 
    } else {
      cf *= (u-d);
    }  
  }
  
  if(h > 0){
    if(abs(c1) > pow(std::numeric_limits<double>::epsilon(), 0.5)){
      cf *= 1.0/c1*(exp(c1 * h) - exp(c1 * l) + exp(-c1*l) - exp(-c1*h));
    } else {
      cf *= 2*(h-l);
    }  
  }
  
  cf -= 1.0;
  
  // Now do the pure vol jump CF
  complex<double> cfvol;
  if(u > 0){
    cfvol = 1.0/(u-d);
    if(abs(c2) > pow(std::numeric_limits<double>::epsilon(), 0.5)){
      cfvol *= 1.0/c2*(exp(c2 * u) - exp(c2 * d)); 
    } else {
      cfvol *= (u-d);
    }  
  } else {
    cfvol = 1.0;
  }
  
  cfvol -= 1.0;
  cf = 0.5 * cf + 0.5 * cfvol;
  return cf;
}

// // [[Rcpp::export]]
// std::complex<double> jumpTransform_1sidedExp(const arma::cx_colvec beta, const Rcpp::List jmpPar){
//   try{
//     double muStock = Rcpp::as<double>(jmpPar["muStock"]); // exponential parameter for negative asset price jumps
//     double muInt = Rcpp::as<double>(jmpPar["muInt"]); // exponential parameter for intensity jump
//     double muVol2 = Rcpp::as<double>(jmpPar["muVol2"]); // exponential parameter for intensity jump / 2nd vol jump
//     double rhoc = Rcpp::as<double>(jmpPar["rhoc"]); // 'correlation' parameter for first vol and stock jump
//     double muVol = Rcpp::as<double>(jmpPar["muVol"]); // exponential parameter for first volatility jump
//     double gammaProp = Rcpp::as<double>(jmpPar["gammaProp"]); // proportionality parameter for balancing jump types against each other
//     
//     // jump in the underlying and the first volatility factor 
//     complex<double> c = muVol * beta(1) + rhoc * muVol * beta(0);
//     
//     complex<double> cf = 1.0 / (1.0 + beta(0) * muStock);
//     cf *= 1.0/(complex<double>(1,0) - c);
//     cf -= complex<double>(1,0);
//     
//     // jump in the intensity factor and other volatility factor
//     cf += gammaProp * (1.0 / (1.0 - beta(2) * muInt) / (1.0 - beta(3) * muVol2) - 1.0);
//     
//     return(cf);
//     
//   } catch( std::exception& __ex__) {
//     forward_exception_to_r(__ex__);
//   } catch(...) {
//     ::Rf_error( "c++ exception (unknown reason)" );
//   }
//   return 0;
// }
// 
// // [[Rcpp::export]]
// std::complex<double> jumpTransform_1sidedExp_2(const arma::cx_colvec beta, const Rcpp::List jmpPar){
//   try{
//     double muStock = Rcpp::as<double>(jmpPar["muStock"]); // exponential parameter for negative asset price co-jumps
//     double muStock2 = Rcpp::as<double>(jmpPar["muStock2"]); // exponential parameter for individual jump in asset price
//     double rhoc = Rcpp::as<double>(jmpPar["rhoc"]); // 'correlation' parameter for first vol and stock jump
//     double muVol = Rcpp::as<double>(jmpPar["muVol"]); // exponential parameter for first volatility jump
//     double gammaProp = Rcpp::as<double>(jmpPar["gammaProp"]); // proportionality parameter for balancing jump types against each other
//     
//     // jump in the underlying and the first volatility factor 
//     complex<double> c = muVol * beta(1) + rhoc * muVol * beta(0);
//     
//     complex<double> cf = 1.0 / (1.0 + beta(0) * muStock);
//     cf *= 1.0/(complex<double>(1,0) - c);
//     cf -= complex<double>(1,0);
//     
//     // jump in the intensity factor and other volatility factor
//     cf += gammaProp * (1.0 / (1.0 + beta(0) * muStock2) - complex<double>(1.0,0));
//     
//     return(cf);
//     
//   } catch( std::exception& __ex__) {
//     forward_exception_to_r(__ex__);
//   } catch(...) {
//     ::Rf_error( "c++ exception (unknown reason)" );
//   }
//   return 0;
// } 