// accurate implementation of the complex expm1 function
#include <complex>
#include "Rcpp.h"

using namespace std;

complex<double> expm1c(complex<double> z)
{
  complex<double> ii(0,1);
  return exp(imag(z)*ii) * expm1(real(z)) + ii * sin(imag(z)) - 2 * sin(imag(z)/2) * sin(imag(z)/2);
}