#include "makePosDef.h" 
#include "RcppArmadillo.h"

using namespace arma;
using namespace std;
using namespace Rcpp;

void makePosDef(mat& semiDefMat, double relEig) {
  try{
      // make positive definite, by first taking an eigenvalue decomposition
      vec eigVal;
      mat eigVec;
      eig_sym(eigVal,eigVec,semiDefMat);
      
      
      // what should the minimal eigenvalue be? Note that the last element of eigVal is the largest. If all eigenvalues are negative, then return the zero matrix.
      if (eigVal.max() < 0) {
        semiDefMat *= 0;
      } else {
        double minEig = eigVal.max() * relEig;
      
        for (int uu=0; uu<eigVal.n_elem - 1; uu++) {
          if (eigVal(uu) < minEig) {
            eigVal(uu) = minEig;
          }
        }
      
        semiDefMat = eigVec * diagmat(eigVal) * trans(eigVec);
      }
  } catch( std::exception& __ex__) {
    forward_exception_to_r(__ex__);
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
}

//' @title Make Positive Definite
//' @description Take a semi-positive definite matrix and use the eigenvalue decomposition to obtain the nearest positive definite matrix
//' @param semiDefMat semi-pos-def matrix, symmetric
//' @param relEig relative range of eigenvalues taken to define the pd matrix
//' @return matrix size of semiDefMat
//' @export
// [[Rcpp::export]]
arma::mat makePositiveDefinite(arma::mat& semiDefMat, double relEig = 1e-6) {
  try{
    // make positive definite, by first taking an eigenvalue decomposition
    arma::vec eigVal;
    arma::mat eigVec;
    arma::eig_sym(eigVal,eigVec,semiDefMat);
    
    
    // what should the minimal eigenvalue be? Note that the last element of eigVal is the largest. If all eigenvalues are negative, then return the zero matrix.
    if (eigVal.max() < 0) {
      semiDefMat *= 0;
    } else {
      double minEig = eigVal.max() * relEig;
      
      for (int uu=0; uu<eigVal.n_elem - 1; uu++) {
        if (eigVal(uu) < minEig) {
          eigVal(uu) = minEig;
        }
      }
      
      semiDefMat = eigVec * diagmat(eigVal) * trans(eigVec);
    }
    return semiDefMat;
  } catch( std::exception& __ex__) {
    forward_exception_to_r(__ex__);
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return 0;
}
