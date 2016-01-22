#include <R_ext/Rdynload.h>
#include <RcppArmadillo.h>
#include <cmath>
#include <complex>
#include "expm1c.h"
#include "jumpTransform.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

static SEXP _RDy_aff_parms;
//static Rcpp::List D;

extern "C" void initmod(void(* odeparms)(int *, double *))
{
  DL_FUNC get_deSolve_gparms;
  SEXP gparms;
  get_deSolve_gparms = R_GetCCallable("deSolve","get_deSolve_gparms");
  gparms = get_deSolve_gparms();
  _RDy_aff_parms = gparms;
}

extern "C" void derivs2 (int *neq, double *t, complex<double> *y, complex<double> *ydot, double *yout, int *ip)
{
  Rcpp::List D(_RDy_aff_parms);
  
  // get the number of factors
  int Nfactors = Rcpp::as<int>(D["N.factors"]);
  
  // initialize complex matrices
  NumericVector h1r = Rcpp::as<NumericVector>(D["H1r"]);
  NumericVector h1i = Rcpp::as<NumericVector>(D["H1i"]);
  cube H1r(h1r.begin(),Nfactors+1,Nfactors+1,Nfactors+1,false);
  cube H1i(h1i.begin(),Nfactors+1,Nfactors+1,Nfactors+1,false);
  cx_cube H1(H1r,H1i);
  
  NumericMatrix k1r = Rcpp::as<NumericMatrix>(D["K1r"]);
  NumericMatrix k1i = Rcpp::as<NumericMatrix>(D["K1i"]);
  mat K1r(k1r.begin(),Nfactors+1,Nfactors+1,false);
  mat K1i(k1i.begin(),Nfactors+1,Nfactors+1,false);
  cx_mat K1(K1r,K1i);
  
  cx_vec K0 = Rcpp::as<arma::cx_vec>(D["K0"]);
  
  // initalize scalars corresponding to jump intensities
  NumericVector L1 = Rcpp::as<NumericVector>(D["l1"]);
  colvec l1 = vec(L1.begin(),Nfactors+1,false);
  
  NumericVector L0 = Rcpp::as<NumericVector>(D["l0"]);
  colvec l0(L0.begin(),1,false);
  
  // now get jump parameters
  double muYc = Rcpp::as<double>(D["muYc"]);
  double sigmaYc = Rcpp::as<double>(D["sigmaYc"]);
  double muSc = Rcpp::as<double>(D["muSc"]);
  double rhoc = Rcpp::as<double>(D["rhoc"]);

  // pick jump transform
  std::string transformName = Rcpp::as<std::string>(D["transformName"]);
  std::complex<double> (*jumpTrFoo)(const arma::cx_colvec&, const double&, const double&, const double&, const double&);
  if(transformName == "expNormJumpTransform"){
    jumpTrFoo = &jumpTransform;
  } else if(transformName == "kouExpJumpTransform"){
    jumpTrFoo = &kouExpTransform;
  } else {
    throw std::range_error("Wrong transform name, check value!");
  }
  
  cx_colvec beta(y,Nfactors+1,false);

  // calculate jump transform
  complex<double> jmpTr = jumpTrFoo(beta,muYc,sigmaYc,muSc,rhoc);
  // complex<double> jmpTr = jumpTransform(beta,muYc,sigmaYc,muSc,rhoc);
  
  cx_colvec yRes = trans(K1) * beta;
  
  for (int i=0;i<Nfactors+1;i++) {
    ydot[i] = yRes[i];
    ydot[i] += l1(i) * jmpTr;
  }
  
  // now add vol of vol
  for (int i=1;i<Nfactors+1;i++) {
    ydot[i] +=  0.5 * ((complex<double>) as_scalar(dot(strans(beta) * H1.slice(i),beta)));
  }
  
  // now do the alpha
  ydot[Nfactors+1] = (complex<double>) dot(K0,beta) + l0(0) * jmpTr;  
}

extern "C" void derivsExt (int *neq, double *t, complex<double> *y, complex<double> *ydot, double *yout, int *ip){
  
  Rcpp::List D(_RDy_aff_parms);
  
  // get the number of factors
  int Nfactors = Rcpp::as<int>(D["N.factors"]);
  
  // initialize complex matrices
  NumericVector h1r = Rcpp::as<NumericVector>(D["H1r"]);
  NumericVector h1i = Rcpp::as<NumericVector>(D["H1i"]);
  cube H1r(h1r.begin(),Nfactors+1,Nfactors+1,Nfactors+1,false);
  cube H1i(h1i.begin(),Nfactors+1,Nfactors+1,Nfactors+1,false);
  cx_cube H1(H1r,H1i);
  
  NumericMatrix k1r = Rcpp::as<NumericMatrix>(D["K1r"]);
  NumericMatrix k1i = Rcpp::as<NumericMatrix>(D["K1i"]);
  mat K1r(k1r.begin(),Nfactors+1,Nfactors+1,false);
  mat K1i(k1i.begin(),Nfactors+1,Nfactors+1,false);
  cx_mat K1(K1r,K1i);
  
  cx_vec K0 = Rcpp::as<arma::cx_vec>(D["K0"]);
  
  // initalize scalars corresponding to jump intensities
  NumericVector L1 = Rcpp::as<NumericVector>(D["l1"]);
  colvec l1 = vec(L1.begin(),Nfactors+1,false);
  
  NumericVector L0 = Rcpp::as<NumericVector>(D["l0"]);
  colvec l0(L0.begin(),1,false);
  
  // now get jump parameters
  double muYc = Rcpp::as<double>(D["muYc"]);
  double sigmaYc = Rcpp::as<double>(D["sigmaYc"]);
  double muSc = Rcpp::as<double>(D["muSc"]);
  double rhoc = Rcpp::as<double>(D["rhoc"]);
  
  // pick jump transform
  std::string transformName = Rcpp::as<std::string>(D["transformName"]);
  std::complex<double> (*jumpTrFoo)(const arma::cx_colvec&, const double&, const double&, const double&, const double&);
  arma::cx_mat (*jumpTrD1Foo)(const arma::cx_colvec&, const double&, const double&, const double&, const double&);
  arma::cx_mat (*jumpTrD2Foo)(const arma::cx_colvec&, const double&, const double&, const double&, const double&);
  arma::cx_mat (*jumpTrD3Foo)(const arma::cx_colvec&, const double&, const double&, const double&, const double&);
  jumpTrFoo = &jumpTransform;
  jumpTrD1Foo = &jumpTransformD1;
  jumpTrD2Foo = &jumpTransformD2;
  jumpTrD3Foo = &jumpTransformD3;
  if(transformName == "expNormJumpTransform"){
  } else if(transformName == "kouExpJumpTransform"){
    jumpTrFoo = &kouExpTransform;
    jumpTrD1Foo = &kouExpTransformD1;
    jumpTrD2Foo = &kouExpTransformD2;
    jumpTrD3Foo = &kouExpTransformD3;
  }
  
  
  // layout of y vector: 5 entries per derivative, total (2+N.factors)*4
  // let kk denote the derivative order
  // stock coeff indexing : kk*(2+Nfactors)
  // constant (alpha) indexing:: kk*(2+Nfactors)+(Nfactors+1)
  // factor coeff indexing: kk*(2+Nfactors)+1 : kk*(2+Nfactors)+Nfactors
  
  cx_colvec yArma(y,4*(Nfactors+2),false);
  cx_colvec beta(Nfactors+1,fill::zeros);
  beta = yArma.subvec(0,Nfactors);
  
  cx_colvec betaPrime(Nfactors+1,fill::zeros);
  betaPrime = yArma.subvec(1*(2+Nfactors), 1*(2+Nfactors)+Nfactors);
  
  cx_colvec betaDblPrime(Nfactors+1,fill::zeros);
  betaDblPrime = yArma.subvec(2*(2+Nfactors),2*(2+Nfactors)+Nfactors);
  
  cx_colvec betaTplPrime(Nfactors+1,fill::zeros);
  betaTplPrime = yArma.subvec(3*(2+Nfactors),3*(2+Nfactors)+Nfactors);
  
  // calculate jump transform
  complex<double> jmpTr;
  jmpTr = jumpTrFoo(beta,muYc,sigmaYc,muSc,rhoc);
  
  // three jump transform derivatives
  cx_mat jmpTrD1 = jumpTrD1Foo(beta,muYc,sigmaYc,muSc,rhoc);
  cx_mat jmpTrD1Tot = jmpTrD1 * betaPrime;
  
  cx_mat jmpTrD2 = jumpTrD2Foo(beta,muYc,sigmaYc,muSc,rhoc);
  cx_mat jmpTrD2Tot = jmpTrD1 * betaDblPrime + betaPrime.st() * jmpTrD2 * betaPrime;
  
  cx_mat jmpTrD3 = jumpTrD3Foo(beta,muYc,sigmaYc,muSc,rhoc);
  
  cx_mat unvecH = jmpTrD3 * betaPrime;
  unvecH.reshape((1+Nfactors),(1+Nfactors));
  
  cx_mat jmpTrD3Tot = jmpTrD1 * betaTplPrime + 3.0 * betaPrime.st() * jmpTrD2 * betaDblPrime + betaPrime.st() * unvecH * betaPrime;
  
  // main ODE: vol driver
  cx_colvec yRes(Nfactors+2);
  yRes.subvec(0,Nfactors) = trans(K1) * beta;
  
  // main ODE: vol-of-vol and jumps
  for (int i=0;i<Nfactors+1;i++) {
    yRes(i) +=  0.5 * ((complex<double>) as_scalar(dot(strans(beta) * H1.slice(i),beta)));
    yRes(i) += l1(i) * jmpTr;
  }
  
  // main ODE: alpha
  yRes(Nfactors+1) = (complex<double>) dot(K0,beta) + l0(0) * jmpTr;
  
  // first deriv ODE
  cx_colvec yResPrime(Nfactors+2);
  yResPrime.subvec(0,Nfactors) = strans(K1) * betaPrime;
  
  // first deriv ODE: vol-of-vol and jumps
  for (int i=0;i<Nfactors+1;i++) {
    yResPrime(i) +=  ((complex<double>) as_scalar(dot(strans(betaPrime) * H1.slice(i),beta)));
    yResPrime(i) += l1(i) * jmpTrD1Tot(0);
  }
  
  // first Deriv: alpha
  yResPrime(Nfactors+1) = (complex<double>) as_scalar(dot(K0, betaPrime));
  yResPrime(Nfactors+1) += l0(0) * jmpTrD1Tot(0);
  
  // second deriv ODE
  cx_colvec yResDblPrime(Nfactors+2);
  yResDblPrime.subvec(0,Nfactors) = strans(K1) * betaDblPrime;
  
  // second deriv ODE: vol-of-vol and jumps
  for (int i=0;i<Nfactors+1;i++) {
    yResDblPrime(i) +=  ((complex<double>) as_scalar(dot(strans(betaDblPrime) * H1.slice(i),beta)));
    yResDblPrime(i) +=  ((complex<double>) as_scalar(dot(strans(betaPrime) * H1.slice(i),betaPrime)));
    yResDblPrime(i) += l1(i) * jmpTrD2Tot(0);
  }
  
  // second Deriv: alpha
  yResDblPrime(Nfactors+1) = (complex<double>) as_scalar(dot(K0, betaDblPrime));
  yResDblPrime(Nfactors+1) += l0(0) * jmpTrD2Tot(0);
  
  // third deriv ODE
  cx_colvec yResTplPrime(Nfactors+2);
  yResTplPrime.subvec(0,Nfactors) = strans(K1) * betaTplPrime;
  // third deriv ODE: vol-of-vol and jumps
  for (int i=0;i<Nfactors+1;i++) {
    yResTplPrime(i) +=  ((complex<double>) as_scalar(dot(strans(betaTplPrime) * H1.slice(i),beta)));
    yResTplPrime(i) +=  3.0*((complex<double>) as_scalar(dot(strans(betaDblPrime) * H1.slice(i),betaPrime)));
    yResTplPrime(i) += l1(i) * jmpTrD3Tot(0);
  }
  
  // third Deriv: alpha
  yResTplPrime(Nfactors+1) = (complex<double>) as_scalar(dot(K0, betaTplPrime));
  yResTplPrime(Nfactors+1) += l0(0) * jmpTrD3Tot(0);
  
  for(int kk=0; kk < Nfactors+2; kk++){
    ydot[kk] = yRes(kk);
    ydot[Nfactors+2+kk] = yResPrime(kk);
    ydot[2*(Nfactors+2)+kk] = yResDblPrime(kk);
    ydot[3*(Nfactors+2)+kk] = yResTplPrime(kk);
  }
}