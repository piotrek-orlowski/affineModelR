#include <RcppArmadillo.h>
#include "../inst/include/affineModelR.h"

using namespace std;
using namespace arma;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List affineSimulateCpp(SEXP TT_, SEXP BB_, SEXP parList_, SEXP dt_, SEXP initVals_, SEXP genPtr_, SEXP retainIndex_){
    try{
      
      // Access the random number generator state
      RNGScope scope;
      
      // convert dt_ to double
      double dt = as<double>(dt_);
      
      // convert variables determining sim size to ints    
      int TT = as<int>(TT_);
      int BB = as<int>(BB_);
      
      // This indicates which observations are retained and which are discarded
      vec retainIndex = as<vec>(retainIndex_);
      
      // Initialise stock and variance factor matrices for return. They have to have one row more than target sim length so that you can initialise and discard nicely.
      int simLength = TT;
      int simWidth = 1;
      // mat sArray = mat(simLength+1, simWidth, fill::zeros);
      // mat vArray = mat(simLength+1,(BB/2)*simWidth, fill::zeros);
      mat sArray = mat(retainIndex.size()+1L, simWidth, fill::zeros);
      mat vArray = mat(retainIndex.size()+1L,(BB/2)*simWidth, fill::zeros);
      mat jumpMarks = mat(retainIndex.size()+1L, 1+(BB/2)*simWidth, fill::zeros);
      // vec jumpValues(1+BB/2,fill::zeros);
      
      // mat sArrayTemp = mat(1, simWidth, fill::zeros);
      // mat vArrayTemp = mat(1, (BB/2.0)*simWidth, fill::zeros);
      
      // Work with the initVals_ list to fill in the first rows of sArray and vArray
      List initVals(initVals_);
      rowvec sInit = as<rowvec>(initVals["S.array"]);
      sArray.row(0) = sInit;
      // sArrayTemp = sInit;
      
      rowvec vInit = as<rowvec>(initVals["V.array"]);
      vArray.row(0) = vInit;
      // vArrayTemp = vInit;
      
      // square-root dt for Euler-scheme purposes
      double dtSqrt = pow(dt,0.5);
      
      // Walk through the parList argument and create vectors which load on stuff in propagation equations
      List parList(parList_);
      
      vec stockDtTerms = vec(1);
      stockDtTerms(0) = as<double>(parList["terms.dt"]);
      colvec stockVdtTerms = as<colvec>(parList["terms.vdt"]);
      colvec stockVdWTerms = as<colvec>(parList["terms.vdW"]);
      colvec stockVdWortTerms = as<colvec>(parList["terms.vdWort"]);
      mat intensity = as<mat>(parList["intensity.terms"]);
      
      mat volDtTerms = as<mat>(parList["vterms.dt"]);
      mat volVdtTerms = as<mat>(parList["vterms.vdt"]);
      mat volVdWTerms = as<mat>(parList["vterms.vdW"]);
      
      // Save number of jumps
      int numJmpTransforms = parList["numJmpTransforms"];
      
      //// Set up state and stock propagation in first-order Euler scheme.
      // initialise holders for past variables and increments
      colvec sPrev = sArray.row(0).t();
      colvec vPrev = vArray.row(0).t();
      colvec vPrevInt = colvec(1+BB/2);
      vPrevInt(0) = 1;
      vPrevInt.rows(1,BB/2) = vPrev;
      colvec dLogS = colvec(sPrev.n_rows);
      colvec dv = colvec(vPrev.n_rows);
      colvec instIntensity = colvec(numJmpTransforms,fill::zeros);
      colvec oneVec = colvec(BB/2,fill::ones);
      colvec vNew = colvec(BB/2);
      
      // Generate vector for jump checks
      // colvec jump = genFoo(parList["jmpPar"]);
      int jmpLength = 1 + vPrev.n_elem;
      
      // Count jumps
      vec numJumps(jmpLength, arma::fill::zeros);
      double cumdt = 0;
      arma::mat jumpSizes(jmpLength, retainIndex.size(), arma::fill::zeros);
      
      // Count iterations so that you know when to retain
      int iterationCounter = 1L;
      
      // for(int ii=1; ii < (TT+1); ii++){
      for(int ii=1; ii < (retainIndex(retainIndex.n_elem-1)+1); ii++){
        
        // count where you are in sample
        cumdt += dt;
        
        // generate random normals
        mat bmGrid = mat(1,BB,fill::randn);
        bmGrid *= dtSqrt;
        
        // set previous values
        // sPrev = sArray.row(ii-1).t();
        // sPrev = sArrayTemp.row(0).t();
        // vPrev = vArray.row(ii-1).t();
        // vPrev = vArrayTemp.row(0).t();
        vPrevInt.rows(1,BB/2) = vPrev;
        // current jump intensity (vector of intensities for all jump types)
        instIntensity = intensity.t() * vPrevInt;
        //// stock increment without jumps
        // constant terms
        dLogS = stockDtTerms;
        // terms multiplying vols that enter the drift
        dLogS += vPrev.t() * stockVdtTerms;
        // multiply by dt
        dLogS = dLogS * dt;
        // Brownian innovations to dLogS, first for BMs correlated with those that drive variance factors
        dLogS += bmGrid.submat(0,0,0,BB/2-1) * (stockVdWTerms % sqrt(vPrev));
        
        // Brownian innovations to dLogS, for BMs uncorrelated with those that drive variance factors
        dLogS += bmGrid.submat(0,BB/2,0,BB-1) * (stockVdWortTerms % sqrt(vPrev));
        
        //// vol increment without jumps
        dv = volDtTerms.t() * oneVec * dt;
        dv += volVdtTerms.t() * vPrev * dt;
        dv += pow(volVdWTerms.t() * vPrev,0.5) % bmGrid.submat(0,0,0,BB/2-1).t();
        
        // Is there a jump? Generate, if yes
        // Check for each factor, jmpProb is a vector wth numJmpTransforms entries
        vec jmpProb = 1.0 - arma::exp(-instIntensity * dt);
        vec jumpCheck(jmpProb.n_elem, arma::fill::randu);
        
        for(int j = 0; j < numJmpTransforms; j++){
          double jumpCheckLoc = jumpCheck(j);
          // Get the jump generator pointer
          Rcpp::List ptrList_ = genPtr_;
          SEXP genPtrLoc = ptrList_[j];
          XPtr<funcPtr> genPtr(genPtrLoc);
          funcPtr genFoo = *genPtr;
          Rcpp::List jmpParLoc = parList["jmpPar"];
          if(jumpCheckLoc < jmpProb(j)){
            numJumps(j)++;
            jumpMarks(iterationCounter, j) = cumdt;
            colvec jump = genFoo(jmpParLoc[j]);
            
            jumpSizes.col(iterationCounter) = jump;
            
            int locNumJumps = jump.n_elem;
            dLogS += jump(0);
            if(locNumJumps > 1){
              dv.subvec(0,locNumJumps-2) += jump.subvec(1,locNumJumps-1);
            }
          }
        }
        // if(jumpCheck(0) < jmpProb){
        //   numJumps++;
        //   jumpMarks(iterationCounter) = cumdt;
        //   jump = genFoo(parList["jmpPar"]);
        //   jumpSizes.col(iterationCounter) += jump;
        //   
        //   dLogS += jump(0);
        //   
        //   dv.subvec(0,jmpLength-2) += jump.subvec(1,jmpLength-1);
        // }
        // add increment to previous value, store
        sPrev = sPrev + dLogS;
        
        vPrev = vPrev + dv;
        if(ii == retainIndex(iterationCounter)){
          sArray.row(iterationCounter) = sPrev;
          if(all(vPrev > 0.0)){
            vArray.row(iterationCounter) = vPrev.t();
          }
          else{
            for(unsigned int kk = 0; kk < vArray.row(iterationCounter).n_cols; kk++){
              vArray(iterationCounter,kk) = max(1e-9,vPrev(kk));
            }
          }
          ++iterationCounter;
        }
        if(!all(vPrev > 0.0)){
          for(unsigned int kk = 0; kk < vArray.row(iterationCounter).n_cols; kk++){
            vPrev(kk) = max(1e-9,vPrev(kk));
          }
        }
      }
      // exponentiate sArray
      sArray = exp(sArray);
      
      // remove last entry from sArray and vArray
      //  they were added because of some strange explosions
      sArray = sArray.head_rows(retainIndex.n_elem);
      vArray = vArray.head_rows(retainIndex.n_elem);
      
      List returnList = List::create(Named("S.array") = sArray, Named("V.array") = vArray, Named("num.jumps") = numJumps, Named("dt") = dt, Named("jump.times") = jumpMarks, Named("jump.sizes") = jumpSizes);
      
      return wrap(returnList);
    } catch( std::exception& __ex__) {
      forward_exception_to_r(__ex__);
    } catch( std::overflow_error& __ex__) {
      forward_exception_to_r(__ex__);
    } catch(...) {
      ::Rf_error( "c++ exception (unknown reason)" );
    } 
    return wrap(0);
  }