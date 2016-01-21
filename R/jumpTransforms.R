#' non-vectorized jumpTransform of the DPS exponential-normal mixture
#' @param uu (N.factors+1) length vector where the jump transform should be evaluated
#' @param jmp a structure describing the jump parameters (muY/sigmaY = mean/sd of log stock jump)
#' @export
#' @references Duffie-Pan-Singleton (2000)
#' @return scalar value corresponding to the jump transform evaluated at the given frequency

expNormJumpTransform <- function(uu,jmp){
  # uu is complex vector of length M
  # the function returns a scalar
  cf = exp(jmp$muYc*uu[1]+.5*jmp$sigmaYc^2*uu[1]^2)/(1-jmp$muSc*uu[2]-jmp$rhoc*jmp$muSc*uu[1])-1
  return(cf)
}

#' @export
kouExpJumpTransform <- function(uu, jmp){
  cf = kouExpTransform(beta = uu, muYc = jmp$muYc, sigmaYc = jmp$sigmaYc, muSc = jmp$muSc, rhoc = jmp$rhoc)
  return(cf)
}