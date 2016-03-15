// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef __affineModelR_RcppExports_h__
#define __affineModelR_RcppExports_h__

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace affineModelR {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("affineModelR", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("affineModelR", "affineModelR_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in affineModelR");
            }
        }
    }

    inline arma::cube affineCFevalCpp(const arma::cube& coeffs, const arma::mat& stateMat, const bool retLog) {
        typedef SEXP(*Ptr_affineCFevalCpp)(SEXP,SEXP,SEXP);
        static Ptr_affineCFevalCpp p_affineCFevalCpp = NULL;
        if (p_affineCFevalCpp == NULL) {
            validateSignature("arma::cube(*affineCFevalCpp)(const arma::cube&,const arma::mat&,const bool)");
            p_affineCFevalCpp = (Ptr_affineCFevalCpp)R_GetCCallable("affineModelR", "affineModelR_affineCFevalCpp");
        }
        RObject __result;
        {
            RNGScope __rngScope;
            __result = p_affineCFevalCpp(Rcpp::wrap(coeffs), Rcpp::wrap(stateMat), Rcpp::wrap(retLog));
        }
        if (__result.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (__result.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(__result).c_str());
        return Rcpp::as<arma::cube >(__result);
    }

    inline arma::cube affineCFderivsEvalCpp(const arma::cube& coeffs, const arma::mat& stateMat) {
        typedef SEXP(*Ptr_affineCFderivsEvalCpp)(SEXP,SEXP);
        static Ptr_affineCFderivsEvalCpp p_affineCFderivsEvalCpp = NULL;
        if (p_affineCFderivsEvalCpp == NULL) {
            validateSignature("arma::cube(*affineCFderivsEvalCpp)(const arma::cube&,const arma::mat&)");
            p_affineCFderivsEvalCpp = (Ptr_affineCFderivsEvalCpp)R_GetCCallable("affineModelR", "affineModelR_affineCFderivsEvalCpp");
        }
        RObject __result;
        {
            RNGScope __rngScope;
            __result = p_affineCFderivsEvalCpp(Rcpp::wrap(coeffs), Rcpp::wrap(stateMat));
        }
        if (__result.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (__result.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(__result).c_str());
        return Rcpp::as<arma::cube >(__result);
    }

}

#endif // __affineModelR_RcppExports_h__
