// Taylor expansion of the jump transform

#include "RcppArmadillo.h"

using namespace arma;
using namespace std;

extern "C" colvec jumpTransformTaylor(const colvec& z, const mat& A, const mat& Ainv, const vec& muY, const vec& sigmaY) {
  try{

  colvec sX = Ainv * z;
  
  colvec jt = z * muY;
  
  colvec scaleVec = 0*sX;
  
  scaleVec += as_scalar(pow(sigmaY,2)+pow(muY,2))/2.0E+0 * sX;
  scaleVec += as_scalar(3*muY*pow(sigmaY,2)+pow(muY,3))/6.0E+0 * pow(sX,2);
  scaleVec += as_scalar(3*pow(sigmaY,4)+6*pow(muY,2)*pow(sigmaY,2)+pow(muY,4))/2.4E+1 * pow(sX,3);
  scaleVec += as_scalar(15*muY*pow(sigmaY,4)+10*pow(muY,3)*pow(sigmaY,2)+pow(muY,5))/1.2E+2 * pow(sX,4);
  scaleVec += as_scalar(15*pow(sigmaY,6)+45*pow(muY,2)*pow(sigmaY,4)+15*pow(muY,4)*pow(sigmaY,2)+pow(muY,6))/7.2E+2 * pow(sX,5);
  scaleVec += as_scalar(105*muY*pow(sigmaY,6)+105*pow(muY,3)*pow(sigmaY,4)+21*pow(muY,5)*pow(sigmaY,2)+pow(muY,7))/5.04E+3 * pow(sX,6);
  scaleVec += as_scalar(105*pow(sigmaY,8)+420*pow(muY,2)*pow(sigmaY,6)+210*pow(muY,4)*pow(sigmaY,4)+28*pow(muY,6)*pow(sigmaY,2)+pow(muY,8))/4.032E+4 * pow(sX,7);
  scaleVec += as_scalar(945*muY*pow(sigmaY,8)+1260*pow(muY,3)*pow(sigmaY,6)+378*pow(muY,5)*pow(sigmaY,4)+36*pow(muY,7)*pow(sigmaY,2)+pow(muY,9))/3.6288E+5 * pow(sX,8);
  scaleVec += as_scalar(945*pow(sigmaY,10)+4725*pow(muY,2)*pow(sigmaY,8)+3150*pow(muY,4)*pow(sigmaY,6)+630*pow(muY,6)*pow(sigmaY,4)+45*pow(muY,8)*pow(sigmaY,2)+pow(muY,10))/3.6288E+6 * pow(sX,9);
  scaleVec += as_scalar(10395*muY*pow(sigmaY,10)+17325*pow(muY,3)*pow(sigmaY,8)+6930*pow(muY,5)*pow(sigmaY,6)+990*pow(muY,7)*pow(sigmaY,4)+55*pow(muY,9)*pow(sigmaY,2)+pow(muY,11))/3.99168E+7 * pow(sX,10);
  scaleVec += as_scalar(10395*pow(sigmaY,12)+62370*pow(muY,2)*pow(sigmaY,10)+51975*pow(muY,4)*pow(sigmaY,8)+13860*pow(muY,6)*pow(sigmaY,6)+1485*pow(muY,8)*pow(sigmaY,4)+66*pow(muY,10)*pow(sigmaY,2)+pow(muY,12))/4.790016E+8 * pow(sX,11);
  
  return(jt + (A * diagmat(scaleVec) * Ainv) * z);
  
  } catch( std::exception& __ex__) {
    forward_exception_to_r(__ex__);
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return 0;
}