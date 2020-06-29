#include <R_ext/Rdynload.h>
#include <RcppArmadillo.h>
#include <cmath>
#include <complex>
#include "../inst/include/affineModelR.h"
#include "../inst/include/expm1c.h"
#include "../inst/include/jumpTransform.h"

#include <typeinfo>

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
  
  // get the number of jumps
  int Njumps = Rcpp::as<int>(D["N.jumps"]);
  
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
  NumericMatrix L1 = Rcpp::as<NumericMatrix>(D["l1"]);
  mat l1_only_vol = mat(L1.begin(), Nfactors, Njumps, false);
  mat l1(Nfactors+1, Njumps, fill::zeros);
  // assign intensities to lower Nfactor rows (top is for stock)
  l1.rows(1,Nfactors) = l1_only_vol;
  
  NumericVector L0 = Rcpp::as<NumericVector>(D["l0"]);
  rowvec l0(L0.begin(), Njumps, false);
  
  // now get jump parameters
  Rcpp::List jmpPar = D["jmpPar"];

  // pick jump transforms
  Rcpp::List ptrList_ = D["jumpTransformPtr"];
  SEXP tmpPtr;
  
  // initialize ODE dependent variable
  cx_colvec beta(y,Nfactors+1,false);

  // calculate jump transforms
  cx_vec jmpTr(Njumps);
  for(int j = 0; j < Njumps; j++){
    tmpPtr = ptrList_[j];
    Rcpp::XPtr<cmpFuncPtr> jumpTrFoo_(tmpPtr);
    cmpFuncPtr jumpTrFoo = *jumpTrFoo_;
    jmpTr(j) = jumpTrFoo(beta, jmpPar[j]); // check
  }
  // complex<double> jmpTr = jumpTrFoo(beta,jmpPar);
  
  cx_colvec yRes = trans(K1) * beta;
  
  for (int i=0;i<Nfactors+1;i++) {
    ydot[i] = yRes[i];
    for(int j=0; j<Njumps;j++){
      ydot[i] += l1(i,j) * jmpTr(j); 
    }
  }
  
  // now add vol of vol
  for (int i=1;i<Nfactors+1;i++) {
    ydot[i] +=  0.5 * ((complex<double>) as_scalar(dot(strans(beta) * H1.slice(i),beta)));
  }
  
  // now do the alpha
  ydot[Nfactors+1] = (complex<double>) dot(K0,beta);  
  for(int j=0; j<Njumps;j++){
    ydot[Nfactors+1] += l0(j) * jmpTr(j); 
  }
}

extern "C" void derivsExt (int *neq, double *t, complex<double> *y, complex<double> *ydot, double *yout, int *ip){
  
  Rcpp::List D(_RDy_aff_parms);
  
  // get the number of factors
  int Nfactors = Rcpp::as<int>(D["N.factors"]);
  
  // get the number of jumps
  int Njumps = Rcpp::as<int>(D["N.jumps"]);
  
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
  NumericMatrix L1 = Rcpp::as<NumericMatrix>(D["l1"]);
  mat l1_only_vol = mat(L1.begin(), Nfactors, Njumps, false);
  mat l1(Nfactors+1, Njumps, fill::zeros);
  // assign intensities to lower Nfactor rows (top is for stock)
  l1.rows(1,Nfactors) = l1_only_vol;
  
  NumericVector L0 = Rcpp::as<NumericVector>(D["l0"]);
  rowvec l0(L0.begin(), Njumps, false);
  
  // now get jump parameters
  Rcpp::List jmpPar = D["jmpPar"];
  
  
  // pick jump transforms
  Rcpp::List ptrList_ = D["jumpTransformPtr"];
  Rcpp::List D1PtrList_ = D["jumpTransformD1Ptr"];
  Rcpp::List D2PtrList_ = D["jumpTransformD2Ptr"];
  Rcpp::List D3PtrList_ = D["jumpTransformD3Ptr"];
  
  // SEXP tmpPtr = tmpPointers["TF"];
  // Rcpp::XPtr<cmpFuncPtr> jumpTrFoo_(tmpPtr);
  // cmpFuncPtr jumpTrFoo = *jumpTrFoo_;
  // 
  // // pick jump transform derivative
  // tmpPtr = tmpPointers["D1"];
  // Rcpp::XPtr<cmpFuncPtrMat> jumpTrD1Foo_(tmpPtr);
  // cmpFuncPtrMat jumpTrD1Foo = *jumpTrD1Foo_;
  // 
  // // pick jump transform second derivative
  // tmpPtr = tmpPointers["D2"];
  // Rcpp::XPtr<cmpFuncPtrMat> jumpTrD2Foo_(tmpPtr);
  // cmpFuncPtrMat jumpTrD2Foo = *jumpTrD2Foo_;
  // 
  // // pick jump transform third derivative
  // tmpPtr = tmpPointers["D3"];
  // Rcpp::XPtr<cmpFuncPtrMat> jumpTrD3Foo_(tmpPtr);
  // cmpFuncPtrMat jumpTrD3Foo = *jumpTrD3Foo_;

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
  SEXP tmpPtr;
  cx_vec jmpTr(Njumps);
  for(int j = 0; j < Njumps; j++){
    tmpPtr = ptrList_[j];
    Rcpp::XPtr<cmpFuncPtr> jumpTrFoo_(tmpPtr);
    cmpFuncPtr jumpTrFoo = *jumpTrFoo_;
    jmpTr(j) = jumpTrFoo(beta, jmpPar[j]); // check
  }
  
  // first jump transform derivative
  cx_mat jmpTrD1(Njumps, Nfactors + 1);
  
  for(int j = 0; j < Njumps; j++){
    tmpPtr = D1PtrList_[j];
    Rcpp::XPtr<cmpFuncPtrMat> jumpTrD1Foo_(tmpPtr);
    cmpFuncPtrMat jumpTrD1Foo = *jumpTrD1Foo_;
    jmpTrD1.row(j) = jumpTrD1Foo(beta, jmpPar[j]);
  }
  // the result is a vector of size Njumps x 1
  cx_mat jmpTrD1Tot = jmpTrD1 * betaPrime;
  
  // second jump transform derivative (Hessian of dim Nf+1 x Nf+1) x Njumps (cube)
  cx_cube jmpTrD2(Nfactors + 1, Nfactors + 1, Njumps);
  cx_mat jmpTrD2Tot(Njumps, 1);
  
  for(int j = 0; j < Njumps; j++){
    tmpPtr = D2PtrList_[j];
    Rcpp::XPtr<cmpFuncPtrMat> jumpTrD2Foo_(tmpPtr);
    cmpFuncPtrMat jumpTrD2Foo = *jumpTrD2Foo_;
    jmpTrD2.slice(j) = jumpTrD2Foo(beta, jmpPar[j]);
    
    jmpTrD2Tot.row(j) = jmpTrD1.row(j) * betaDblPrime + betaPrime.st() * jmpTrD2.slice(j) * betaPrime;
  }
  
  // cx_mat jmpTrD2 = jumpTrD2Foo(beta,jmpPar);
  // cx_mat jmpTrD2Tot = jmpTrD1 * betaDblPrime + betaPrime.st() * jmpTrD2 * betaPrime;
  
  // Third jump transform derivative (gradient of vectorized Hessian) (Nf+1)^2 x (Nf+1) x Njumps (cube)
  cx_cube jmpTrD3(pow(Nfactors + 1, 2), Nfactors + 1, Njumps);
  cx_mat jmpTrD3Tot(Njumps, 1);
  
  for(int j = 0; j < Njumps; j++){
    tmpPtr = D3PtrList_[j];
    Rcpp::XPtr<cmpFuncPtrMat> jumpTrD3Foo_(tmpPtr);
    cmpFuncPtrMat jumpTrD3Foo = *jumpTrD3Foo_;
    jmpTrD3.slice(j) = jumpTrD3Foo(beta, jmpPar[j]);
    
    cx_mat unvecH = jmpTrD3.slice(j) * betaPrime;
    unvecH.reshape((1+Nfactors),(1+Nfactors));
    
    jmpTrD3Tot.row(j) = jmpTrD1.row(j) * betaTplPrime + 3.0 * betaPrime.st() * jmpTrD2.slice(j) * betaDblPrime + betaPrime.st() * unvecH * betaPrime;
  }
  
  // cx_mat jmpTrD3 = jumpTrD3Foo(beta,jmpPar);
  // 
  // cx_mat unvecH = jmpTrD3 * betaPrime;
  // unvecH.reshape((1+Nfactors),(1+Nfactors));
  // 
  // cx_mat jmpTrD3Tot = jmpTrD1 * betaTplPrime + 3.0 * betaPrime.st() * jmpTrD2 * betaDblPrime + betaPrime.st() * unvecH * betaPrime;
  
  // main ODE: vol driver
  cx_colvec yRes(Nfactors+2);
  yRes.subvec(0,Nfactors) = trans(K1) * beta;
  
  // main ODE: vol-of-vol and jumps
  for (int i=0;i<Nfactors+1;i++) {
    yRes(i) +=  0.5 * ((complex<double>) as_scalar(dot(strans(beta) * H1.slice(i),beta)));
    for(int j=0;j<Njumps;j++){
      yRes(i) += l1(i,j) * jmpTr(j); 
    }
  }
  
  // main ODE: alpha
  yRes(Nfactors+1) = (complex<double>) dot(K0,beta);
  for(int j=0;j<Njumps;j++){
    yRes(Nfactors+1) += l0(j) * jmpTr(j); 
  }
  
  // first deriv ODE
  cx_colvec yResPrime(Nfactors+2);
  yResPrime.subvec(0,Nfactors) = strans(K1) * betaPrime;
  
  // first deriv ODE: vol-of-vol and jumps
  for (int i=0;i<Nfactors+1;i++) {
    yResPrime(i) +=  ((complex<double>) as_scalar(dot(strans(betaPrime) * H1.slice(i),beta)));
    for(int j=0;j<Njumps;j++){
      yResPrime(i) += l1(i,j) * jmpTrD1Tot(j); 
    }
  }
  
  // first Deriv: alpha
  yResPrime(Nfactors+1) = (complex<double>) as_scalar(dot(K0, betaPrime));
  // yResPrime(Nfactors+1) += l0(0) * jmpTrD1Tot(0);
  for(int j=0;j<Njumps;j++){
    yResPrime(Nfactors+1) += l0(j) * jmpTrD1Tot(j); 
  }
  
  // second deriv ODE
  cx_colvec yResDblPrime(Nfactors+2);
  yResDblPrime.subvec(0,Nfactors) = strans(K1) * betaDblPrime;
  
  // second deriv ODE: vol-of-vol and jumps
  for (int i=0;i<Nfactors+1;i++) {
    yResDblPrime(i) +=  ((complex<double>) as_scalar(dot(strans(betaDblPrime) * H1.slice(i),beta)));
    yResDblPrime(i) +=  ((complex<double>) as_scalar(dot(strans(betaPrime) * H1.slice(i),betaPrime)));
    for(int j=0;j<Njumps;j++){
      yResDblPrime(i) += l1(i,j) * jmpTrD2Tot(j); 
    }
  }
  
  // second Deriv: alpha
  yResDblPrime(Nfactors+1) = (complex<double>) as_scalar(dot(K0, betaDblPrime));
  for(int j=0;j<Njumps;j++){
    yResDblPrime(Nfactors+1) += l0(j) * jmpTrD2Tot(j); 
  }
  
  // third deriv ODE
  cx_colvec yResTplPrime(Nfactors+2);
  yResTplPrime.subvec(0,Nfactors) = strans(K1) * betaTplPrime;
  // third deriv ODE: vol-of-vol and jumps
  for (int i=0;i<Nfactors+1;i++) {
    yResTplPrime(i) +=  ((complex<double>) as_scalar(dot(strans(betaTplPrime) * H1.slice(i),beta)));
    yResTplPrime(i) +=  3.0*((complex<double>) as_scalar(dot(strans(betaDblPrime) * H1.slice(i),betaPrime)));
    for(int j=0;j<Njumps;j++){
      yResTplPrime(i) += l1(i,j) * jmpTrD3Tot(j); 
    }
  }
  
  // third Deriv: alpha
  yResTplPrime(Nfactors+1) = (complex<double>) as_scalar(dot(K0, betaTplPrime));
  for(int j=0;j<Njumps;j++){
    yResTplPrime(Nfactors+1) += l0(j) * jmpTrD3Tot(j); 
  }
  
  for(int kk=0; kk < Nfactors+2; kk++){
    ydot[kk] = yRes(kk);
    ydot[Nfactors+2+kk] = yResPrime(kk);
    ydot[2*(Nfactors+2)+kk] = yResDblPrime(kk);
    ydot[3*(Nfactors+2)+kk] = yResTplPrime(kk);
  }
}