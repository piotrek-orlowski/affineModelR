#include <RcppArmadillo.h>
#include "../inst/include/affineModelR.h"

using namespace std;
using namespace arma;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List affineSimulateCpp(SEXP TT_, SEXP BB_, SEXP parList_, SEXP dt_, SEXP initVals_, SEXP genPtr_){
    try{
      
      // Access the random number generator state
      RNGScope scope;
      
      // Get the jump generator pointer
      XPtr<funcPtr> genPtr(genPtr_);
      funcPtr genFoo = *genPtr;
      
      // convert dt_ to double
      double dt = as<double>(dt_);
      
      // convert variables determining sim size to ints    
      int TT = as<int>(TT_);
      int BB = as<int>(BB_);
      
      // initialise normal random draw matrix
      mat bmGrid = mat(TT,BB,fill::randn);
      
      // Initialise stock and variance factor matrices for return. They have to have one row more than target sim length so that you can initialise and discard nicely.
      int simLength = TT;
      int simWidth = 1;
      mat sArray = mat(simLength+1, simWidth, fill::zeros);
      mat vArray = mat(simLength+1,(BB/2)*simWidth, fill::zeros);
      vec jumpMarks = vec(simLength+1, fill::zeros);
      vec jumpValues(1+BB/2,fill::zeros);
      
      // Work with the initVals_ list to fill in the first rows of sArray and vArray
      List initVals(initVals_);
      rowvec sInit = as<rowvec>(initVals["S.array"]);
      sArray.row(0) = sInit;
      
      rowvec vInit = as<rowvec>(initVals["V.array"]);
      vArray.row(0) = vInit;
      
      // square-root dt for Euler-scheme purposes
      double dtSqrt = pow(dt,0.5);
      //    Rcout << "dt, dtsqrt" << dt << " " << dtSqrt << "\n";
      // scale Normal increments by sqrt(dt)
      bmGrid *= dtSqrt;
      // Rcout << bmGrid;
      
      // Walk through the parList argument and create vectors which load on stuff in propagation equations
      List parList(parList_);
      
      vec stockDtTerms = vec(1);
      stockDtTerms(0) = as<double>(parList["terms.dt"]);
      colvec stockVdtTerms = as<colvec>(parList["terms.vdt"]);
      colvec stockVdWTerms = as<colvec>(parList["terms.vdW"]);
      colvec stockVdWortTerms = as<colvec>(parList["terms.vdWort"]);
      colvec intensity = as<colvec>(parList["intensity.terms"]);
      
      mat volDtTerms = as<mat>(parList["vterms.dt"]);
      mat volVdtTerms = as<mat>(parList["vterms.vdt"]);
      mat volVdWTerms = as<mat>(parList["vterms.vdW"]);
      
      //// Set up state and stock propagation in first-order Euler scheme.
      // initialise holders for past variables and increments
      colvec sPrev = sArray.row(0).t();
      colvec vPrev = vArray.row(0).t();
      colvec vPrevInt = colvec(1+BB/2);
      vPrevInt(0) = 1;
      vPrevInt.rows(1,BB/2) = vPrev;
      colvec dLogS = colvec(sPrev.n_rows);
      colvec dv = colvec(vPrev.n_rows);
      colvec instIntensity = colvec(1,fill::zeros);
      colvec oneVec = colvec(BB/2,fill::ones);
      colvec vNew = colvec(BB/2);
      
      // Generate vector for jump checks
      NumericVector jumpCheckRcpp = runif(TT,0,1);
      colvec jumpCheck(jumpCheckRcpp);
      colvec jmpPar = as<colvec>(parList["jmpPar"]);
      colvec jump = genFoo(jmpPar);
      int jmpLength = jump.n_elem;
      
      double jmpProb = 0.0;
      
      // Count jumps
      int numJumps = 0;
      
      for(int ii=1; ii < (TT+1); ii++){
        // set previous values
        sPrev = sArray.row(ii-1).t();
        vPrev = vArray.row(ii-1).t();
        vPrevInt.rows(1,BB/2) = vPrev;
        // current jump intensity
        instIntensity = intensity.t() * vPrevInt;
        //// stock increment without jumps
        // constant terms
        dLogS = stockDtTerms;
        // terms multiplying vols that enter the drift
        dLogS += vPrev.t() * stockVdtTerms;
        // multiply by dt
        dLogS = dLogS * dt;
        // Brownian innovations to dLogS, first for BMs correlated with those that drive variance factors
        dLogS += bmGrid.submat(ii-1,0,ii-1,BB/2-1) * (stockVdWTerms % sqrt(vPrev));
        //Rcout << "Line 216 dLogS: " << dLogS << "\n";
        // Brownian innovations to dLogS, for BMs uncorrelated with those that drive variance factors
        dLogS += bmGrid.submat(ii-1,BB/2,ii-1,BB-1) * (stockVdWortTerms % sqrt(vPrev));
        
        //// vol increment without jumps
        dv = volDtTerms.t() * oneVec * dt;
        dv += volVdtTerms.t() * vPrev * dt;
        dv += pow(volVdWTerms.t() * vPrev,0.5) % bmGrid.submat(ii-1,0,ii-1,BB/2-1).t();
        
        // Is there a jump? Generate, if yes
        jmpProb = 1.0 - exp(-instIntensity(0) * dt);
        if(jumpCheck(ii-1) < jmpProb){
          numJumps++;
          jumpMarks(ii) = 1.0;
          jump = genFoo(jmpPar);
          //        Rcout << "at ii " << ii << " jmp prob is " << jmpProb << ", stock jump is " << stockJump << ", vol jump is " << volJump <<"\n";
          dLogS += jump(0);
          // Rcout << "jmpLength: \n" << jmpLength << "\n";
          // jump.print("jump itself");
          // Rcout << "dv before a jump: \n" << dv << "\n";
          dv.subvec(0,jmpLength-2) += jump.subvec(1,jmpLength-1);
          // Rcout << "dv with a jump: \n" << dv << "\n";
        }
        // add increment to previous value, store
        sArray.row(ii) = sPrev + dLogS;
        // check whether we have negative vols etc.
        vNew = vPrev + dv;
        if(all(vNew > 0.0)){
          vArray.row(ii) = vNew.t();
        }
        else{
          for(int kk = 0; kk < vArray.row(ii).n_cols; kk++){
            //          vArray(ii,kk) = max(1e-9,vArray(ii,kk));
            vArray(ii,kk) = max(1e-9,vNew(kk));
          }
        }
        //      Rcout << "line 241 vArray row \n" << vArray.row(ii) << "\n";
      }
      //    Rcout << "I'm here, line 246\n";
      //    Rcout << "sArray \n" << sArray << "\n vArray \n" << vArray << "\n"; 
      // exponentiate sArray
      sArray = exp(sArray);
      List returnList = List::create(Named("S.array") = sArray, Named("V.array") = vArray, Named("num.jumps") = numJumps, Named("dt") = dt, Named("jump.times") = jumpMarks);
      
      return wrap(returnList);
      //    return wrap(1);
    } catch( std::exception& __ex__) {
      forward_exception_to_r(__ex__);
    } catch( std::overflow_error& __ex__) {
      forward_exception_to_r(__ex__);
    } catch(...) {
      ::Rf_error( "c++ exception (unknown reason)" );
    } 
    return wrap(0);
  }