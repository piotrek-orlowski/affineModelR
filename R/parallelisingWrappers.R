#' @title Parallelisation with external pointers
#' @description Parallelising calls to R functions in this package while using the pointer-passing interface requires some care: the \code{\link{getPointerToJumpTransform}} and \code{\link{getPointerToGenerator}} functions have to be called \emph{on the cluster}; a call such as \code{parLapply(cl = myCluster, X = uMat.list, fun = affineCF, params.Q = params.Q, params.P = params.P, t.vec = 0.15, v.0 = matrix(0.025,1,1), N.factors = 1, mod.type ='standard', jumpTransform = getPointerToJumpTransform('kouExpJumpTransform'))} will fail because of invalid pointers. This interface delays evaluation of the pointer functions by wrapping the call appropriately.
#' @param foo function to be wrapped
#' @param ptrFoo function that finds pointer, takes a single argument
#' @param ptrArg argument to ptrFoo
#' @param argName character, name of the pointer argument in foo; for example \code{'jumpTransform'} if you want to call \code{\link{affineCF}}
#' @param ptrElem character, if \code{ptrFoo} returns a named list, get pointer in field \code{ptrFoo(ptrArg)$ptrElem}.
#' @details Further functions are ready parallel versions of the CF callers.
#' @return Returns a function that which should be called with the same arguments as \code{foo} except for the pointer arguments.
#' @export

parallelWrapper <- function(foo, ptrFoo, ptrArg, argName, ptrElem = NULL, ...){
  if(is.null(ptrElem)){
    eval(parse(text = paste0('resFoo <- function(...){
        res <- foo(..., ' , argName ,' = ptrFoo(ptrArg))
        return(res)
      }
      return(resFoo)')))
  } else {
    eval(parse(text = paste0('resFoo <- function(...){
        res <- foo(..., ' , argName ,' = ptrFoo(ptrArg)$',ptrElem,')
                             return(res)
        }
        return(resFoo)')))
  }
}