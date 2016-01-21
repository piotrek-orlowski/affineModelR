# Tests the ODE solver correctness given Q parameters
# 
# Author: Andras Sali, Piotr Orlowski
###############################################################################

###### check if black-scholes cf values are correctly returned ----
library(affineModelR)
library(testthat)
N.factors <- 1

params <- list()
params[["1"]] <- list(phi=1,eta=.04,lmb=1e-6,kpp=100,rho=0)
params$jmp <- list(muSc = 1e-8, muYc = 0, sigmaYc=1e-6,lprop=0,lvec=0,rhoc=0)

uVec <- 1:10

mkt <- list(p=1,q=-.5,r=.2,t=.15)

bsCF <- jumpDiffusionODEs(cbind(matrix(uVec,ncol=1)*1i,0),params=params,mkt = mkt, N.factors=N.factors)

bsCF.true <- exp((mkt$r-mkt$q-1/2*params[["1"]]$eta)*mkt$t*uVec*1i - 1/2*params[["1"]]$eta*uVec^2*mkt$t)
names(bsCF.true) <- paste("u",uVec,sep=".")

test_that("ODE solver works for the Black-Scholes model", {
			expect_that(dim(bsCF), equals(c(length(uVec),1,1+N.factors)))
			expect_that(exp(bsCF[,,"a"]+bsCF[,,"b1"]*params[["1"]]$eta), equals(bsCF.true))
		})

#### now make sure that irrelevant factor is indeed irrelevant
params <- list()
params[["1"]] <- list(phi=1,eta=.04,lmb=1e-6,kpp=100,rho=0)
params[["2"]] <- list(phi=0,eta=.04,lmb=1e-1,kpp=1,rho=0)
params$jmp <- list(muSc = 1e-8, muYc = 0, sigmaYc=1e-6,lprop=c(0,0),lvec=0,rhoc=0)

bsCF.2 <- jumpDiffusionODEs(cbind(matrix(uVec,ncol=1)*1i,0,0),params=params,mkt = mkt, N.factors=2)

test_that("ODE solver works for the Black-Scholes model", {
  expect_that(exp(bsCF.2[,,"a"]+bsCF.2[,,"b1"]*params[["1"]]$eta + bsCF.2[,,"b2"]*params[["2"]]$eta), equals(bsCF.true))
})

#### now check if stationary distribution of volatility factor is correct ----
params <- list()
params[["1"]] <- list(phi=1,eta=.04,lmb=.2,kpp=5,rho=0)
params$jmp <- list(muSc = 1e-8, muYc = 0, sigmaYc=1e-6,lprop=c(0),lvec=0,rhoc=0)

mkt.long <- mkt
mkt.long$t <- 20

sv.cf <- jumpDiffusionODEs(matrix(c(rep(0,10),1i*(1:10)),10,2),params=params,mkt = mkt.long, N.factors=1)

cir.shape <- 2*params[["1"]]$kpp*params[["1"]]$eta/params[["1"]]$lmb^2
cir.scale <- params[["1"]]$lmb^2/(2*params[["1"]]$kpp) 

library('prob')

stationary.cf <- cfgamma(1:10,shape=cir.shape,scale=cir.scale)
names(stationary.cf) <- paste("u",uVec,sep=".")

test_that("ODE solver works for the stationary distribution of the first factor", {
			expect_that(exp(sv.cf[,,"a"]), equals(stationary.cf))
		})

#### Now check factor equivalence
dT <- 1/12
mkt <- list(t=dT,p=1,q=0,r=0)

u <- 0

params <- list()
params[["1"]] <- list(phi=1,eta=.04,lmb=1e-1,kpp=10,rho=0)
params[["2"]] <- list(phi=1,eta=.02,lmb=1e-1,kpp=1,rho=-.6)
params$jmp <- list(muSc = 1e-1, muYc = -.1, sigmaYc=.2,lprop=c(0,0),lvec=0,rhoc=-.6)

ode.sol.1 <- matrix((jumpDiffusionODEs(matrix(c(u,1,0),ncol=3),params,mkt=mkt,rtol=1e-15,atol=1e-40,N.factors=2)),nrow=1)

params <- list()
params[["2"]] <- list(phi=1,eta=.04,lmb=1e-1,kpp=10,rho=0)
params[["1"]] <- list(phi=1,eta=.02,lmb=1e-1,kpp=1,rho=-.6)
params$jmp <- list(muSc = 1e-1, muYc = -.1, sigmaYc=.2,lprop=c(0,0),lvec=0,rhoc=-.6)

ode.sol.2 <- matrix((jumpDiffusionODEs(matrix(c(u,0,1),ncol=3),params,mkt=mkt,rtol=1e-15,atol=1e-40,N.factors=2)),nrow=1)

test_that("Factor swap equivalence",{expect_equal(abs(ode.sol.2-ode.sol.1[,c(2,1,3)]),matrix(0,nrow=1,ncol=3),tol=1e-13)})
