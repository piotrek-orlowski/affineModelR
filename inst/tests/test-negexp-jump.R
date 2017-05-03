library(affineModelR)
library(numDeriv)

par.list <- list()
par.list$P <- list()
par.list$P$`1` <- list(erp=0, erp0=0, phi = 0.1, rho = -0.8, kpp = 5, lmb = 1, eta = 1)
par.list$P$`2` <- list(erp=0, phi = 0.0, rho = 0.0, kpp = 3, lmb = 1e-4, eta = 1)
par.list$P$`3` <- list(erp=0, erp0=0, phi = 0.1, rho = 0.0, kpp = 12, lmb = 0.1, eta = 1)
par.list$P$jmp <- list(lvec = 0, lprop = c(0,5,0), gammaProp = 2, muInt = 0.05, muStock = 0.03, muVol = 0.12,  muVol2 = 0.12, rhoc = -0.05)

# par.list.sim <- par.list
# par.list.sim$P$jmp$lprop <- (1 + par.list$P$jmp$gammaProp) * par.list.sim$P$jmp$lprop
# par.list.sim$Q <- par.list.sim$P

# par.list$P$jmp$lvec <- par.list$P$jmp$lvec / (1 + par.list$P$jmp$gammaProp)
# par.list$P$jmp$lprop <- par.list$P$jmp$lprop / (1 + par.list$P$jmp$gammaProp)
par.list$Q <- par.list$P

par.list.kou <- par.list
par.list.kou$Q$jmp <- list(lvec = 0, lprop = c(0,0,0), muYc = -0.0, sigmaYc = 0.0, muSc = 1/0.12, rhoc = -0.0)
par.list.kou$P <- par.list.kou$Q

tf_str <- 'oneSidedExponential'

## sim
mkt <- data.frame(t=10/252,r=0.0,q=0,p=1)
v.0 <- rbind(c(1,1,1),c(3,3,3))

v.0 <- v.0[1,,drop=F]

sim.paths <- affineSimulate(paramsList = par.list, N.factors = 3, t.days = mkt$t * 252, t.freq = 1, freq.subdiv = 86400, rng.seed = 1, jumpGeneratorPtr = getPointerToGenerator(tf_str), jumpTransformPtr = getPointerToJumpTransform(tf_str)$TF, mod.type = 'standard', nrepl = 500, init.vals = list(S.array = 0, V.array = v.0[1,], day.offset = 0), rf.rate = 0)

sim.paths <- affineSimulate(paramsList = par.list, N.factors = 3, t.days = 252, t.freq = 1/78, freq.subdiv = 300, rng.seed = 1, jumpGeneratorPtr = getPointerToGenerator(tf_str), jumpTransformPtr = getPointerToJumpTransform(tf_str)$TF, mod.type = 'standard', nrepl = 1, init.vals = list(S.array = 0, V.array = v.0[1,], day.offset = 0), rf.rate = 0)


sim.paths.kou <- affineSimulate(paramsList = par.list.kou, N.factors = 3, t.days = mkt$t * 252, t.freq = 1, freq.subdiv = 86400, rng.seed = 1, jumpGeneratorPtr = getPointerToGenerator('kouExpJumpTransform'), jumpTransformPtr = getPointerToJumpTransform('kouExpJumpTransform')$TF, mod.type = 'standard', nrepl = 500, init.vals = list(S.array = 0, V.array = v.0[1,], day.offset = 0))


mean(sapply(sim.paths$sim.arrays, function(x) log(tail(x$S.array,1))))
var(sapply(sim.paths$sim.arrays, function(x) log(tail(x$S.array,1))))

mean(sapply(sim.paths.kou$sim.arrays, function(x) log(tail(x$S.array,1))))
var(sapply(sim.paths.kou$sim.arrays, function(x) log(tail(x$S.array,1))))


grad(func = function(x) affineCF(u = matrix(c(x,rep(0,3)),nrow=1), params.Q = par.list$Q, params.P = par.list$P, t.vec = NULL, v.0 = v.0[1,,drop=F], jumpTransform = getPointerToJumpTransform(tf_str)$TF, N.factors = 3, CGF = FALSE, mkt = mkt, atol = 1e-28, rtol = 1e-13)[1], x = 0, method = "complex")

grad(func = function(x) affineCF(u = matrix(c(x,rep(0,3)),nrow=1), params.Q = par.list.kou$Q, params.P = NULL, t.vec = NULL, v.0 = v.0[1,,drop=F], jumpTransform = getPointerToJumpTransform('kouExpJumpTransform')$TF, N.factors = 3, CGF = FALSE, mkt = mkt, atol = 1e-28, rtol = 1e-13)[1], x = 0, method = "complex")

hessian(func = function(x) affineCF(u = matrix(c(x,rep(0,3)),nrow=1), params.Q = par.list$Q, params.P = NULL, t.vec = NULL, v.0 = v.0[1,,drop=F], jumpTransform = getPointerToJumpTransform(tf_str)$TF, N.factors = 3, CGF = FALSE, mkt = mkt, atol = 1e-22, rtol = 1e-13)[1], x = 0, method = "complex") - grad(func = function(x) affineCF(u = matrix(c(x,rep(0,3)),nrow=1), params.Q = par.list$Q, params.P = NULL, t.vec = NULL, v.0 = v.0[1,,drop=F], jumpTransform = getPointerToJumpTransform(tf_str)$TF, N.factors = 3, CGF = FALSE, mkt = mkt, atol = 1e-22, rtol = 1e-13)[1], x = 0, method = "complex")^2

hessian(func = function(x) affineCF(u = matrix(c(x,rep(0,3)),nrow=1), params.Q = par.list.kou$Q, params.P = NULL, t.vec = NULL, v.0 = v.0[1,,drop=F], jumpTransform = getPointerToJumpTransform('kouExpJumpTransform')$TF, N.factors = 3, CGF = FALSE, mkt = mkt, atol = 1e-28, rtol = 1e-13)[1], x = 0, method = "complex") - grad(func = function(x) affineCF(u = matrix(c(x,rep(0,3)),nrow=1), params.Q = par.list.kou$Q, params.P = NULL, t.vec = NULL, v.0 = v.0[1,,drop=F], jumpTransform = getPointerToJumpTransform('kouExpJumpTransform')$TF, N.factors = 3, CGF = FALSE, mkt = mkt, atol = 1e-28, rtol = 1e-13)[1], x = 0, method = "complex")^2


# Plot
layout(t(c(1:4)))
plot(sim.paths$times, sim.paths$sim.arrays[[1]]$S.array, type = 'l', xlab = 'day', ylab = 'Stock', main = 'Heston model \n stock price simulation')
plot(sim.paths$times, sim.paths$sim.arrays[[1]]$V.array[,1], type = 'l', xlab = 'day', ylab = 'Vol', main = 'Heston model \n volatility simulation')
plot(sim.paths$times, sim.paths$sim.arrays[[1]]$V.array[,2], type = 'l', xlab = 'day', ylab = 'Vol', main = 'Heston model \n volatility simulation')
plot(sim.paths$times, sim.paths$sim.arrays[[1]]$V.array[,3], type = 'l', xlab = 'day', ylab = 'Vol', main = 'Heston model \n volatility simulation')

u.ret <- cbind(seq(-2,2,by=0.01),0,0,0) # The first argument is for the stock return. Pad with zeros for volatility.

t.vec <- 1 # plotting horizon is 1 Year

mkt <- data.frame(t=t.vec,r=0.0,q=0,p=1/252)

v.0 <- rbind(c(1,1,1),c(3,3,3))

# We will plot the risk-neutral MGF, thus params.P = NULL.
cf.ret <- affineCF(u = u.ret, params.Q = par.list$Q, params.P = NULL, t.vec = NULL, v.0 = v.0, jumpTransform = getPointerToJumpTransform(tf_str)$TF, N.factors = 3, CGF = FALSE, mkt = mkt)

cf.ret.kou <- affineCF(u = u.ret, params.Q = par.list.kou$Q, params.P = NULL, t.vec = NULL, v.0 = v.0, jumpTransform = getPointerToJumpTransform('kouExpJumpTransform')$TF, N.factors = 3, CGF = FALSE, mkt = mkt)

# Take the real part
cf.ret <- Re(cf.ret)
layout(1)
plot(u.ret[,1], cf.ret[,,1], type = 'l', ylim = range(cf.ret), col = "darkorange", lwd = 2.5, xlab = 'ret arg', ylab = 'ret CF 1Y', main = "return CF -- Heston model")

