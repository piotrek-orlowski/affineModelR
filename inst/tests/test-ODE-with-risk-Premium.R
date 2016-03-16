# This script test simulator and ODE solver consistency for risk-premium returns.
# print("Checking P ODE solver and simulator consistency...")

bsStruct <- namesToStruct()

svFast <- bsStruct$svFast
svFast$lmb <- .2
svFast$kpp <- 5
svFast$rho <- -.1

svSlow <- bsStruct$svSlow
svSlow$lmb <- .1
svSlow$kpp <- 2
svSlow$rho <- -.6
svSlow$eta <- .005

jmp <- bsStruct$jmp
jmp$muY <- -.05
jmp$sigmaY <- .05
jmp$lvec <- c(1,0,0)

params.P <- params.Q <- bsStruct
params.P$svFast <- params.Q$svFast <- svFast
params.P$svSlow <- params.Q$svSlow <- svSlow
params.P$jmp <- params.Q$jmp <- jmp

# now add some risk-premia
params.Q$jmp$muY <- -.1
params.Q$jmp$sigmaY <- .1

params.Q$svFast$kpp <- 3
params.Q$svFast$eta <- .015

params.Q$svSlow$kpp <- 4
params.Q$svSlow$eta <- .01

# set up the parameters controlling the simulation (number of runs, etc..)
N.runs <- 1
uVec <- seq(-3,3)
th.lt <- var.emp.lt <- emp.lt <- 0 * uVec
N.sim <- 800
frac.time <- 1/2

mkt <- data.frame(p=1, r=0, q=0, t= 250/252/frac.time)
# calculate one year conditional lt
th.cond <- twoFactorJumpODEsSolveP(cbind(matrix(uVec,length(uVec),1),0,0), params.P, params.Q,mkt=mkt)

mkt.long <- data.frame(p=1, r=0, q=0, t= 20)

# calculate theoretical moment conditions
for (u in 1:length(uVec)) {
  state.curr <- matrix(c(0,th.cond[u,1,c("b1","b2")]),1,3)
  th.lt[u] <- exp(th.cond[u,1,"a"])*exp(twoFactorJumpODEsSolveP(state.curr, params.P, params.Q,mkt=mkt.long)[,,"a"])
}

# calculate empirical moment conditions, displaying mean + confidence band at each iteration
for (nn in 1:N.runs) {
  #print(paste("Doing run:",nn))
  test <- simulator2fGenerateProcess(params.P,params.Q,t.days=250/frac.time , freq.subdiv=32,t.freq=1,num.paths=N.sim,rng.seed=as.numeric(Sys.time()))
  
  for (u in uVec) {
    emp.lt[which(u==uVec)] <- emp.lt[which(u==uVec)]  + mean((test$S.array[250/frac.time,1:N.sim]/test$S.array[1,1:N.sim])^u)
    var.emp.lt[which(u==uVec)] <- var.emp.lt[which(u==uVec)]  + mean((test$S.array[250/frac.time,1:N.sim]/test$S.array[1,1:N.sim])^(2*u))
  }
  emp.mean <- emp.lt/nn
  var.emp.mean <- var.emp.lt/nn - emp.mean^2
  emp.mean.lower <- emp.mean-3*sqrt(var.emp.mean)/sqrt(nn*N.sim)
  emp.mean.upper <- emp.mean+3*sqrt(var.emp.mean)/sqrt(nn*N.sim)
  #print(rbind(uVec,th.lt,emp.mean,emp.mean.lower,emp.mean.upper))
}
expect_that(all(emp.mean.lower <= th.lt), is_true())
expect_that(all(emp.mean.upper >= th.lt), is_true())
