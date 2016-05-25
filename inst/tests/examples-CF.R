### Examples of calculations possible in the affineModelR package.

# ---- LOAD LIBRARY ----

library(affineModelR)

# ---- LOAD MODEL PARAMETERS ----

# This is Heston's (93) model, reparametrised.
data("heston-parameters")

# ---- PLOT MOMENT-GENERATING FUNCTION ----

# The Heston model has a single stochastic volatility factor. The joint characteristic function of the stock return and vol factor takes two arguments. Here we will plot the marginal MGFs.

u.ret <- cbind(seq(-2,2,by=0.01),0) # The first argument is for the stock return. Pad with zeros for volatility.

t.vec <- 1 # plotting horizon is 1 Year

v.0 <- matrix(c(1,3),nrow = 2,ncol = 1) # We will plot the MGF at factor values 1 (the long-run mean) and 3 (a very high value)

# We will plot the risk-neutral MGF, thus params.P = NULL.
cf.ret <- affineCF(u = u.ret, params.Q = parListHeston$Q, params.P = NULL, t.vec = t.vec, v.0 = v.0, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.factors = 1, CGF = FALSE)

# Take the real part
cf.ret <- Re(cf.ret)

plot(u.ret[,1], cf.ret[,,1], type = 'l', ylim = range(cf.ret), col = "darkorange", lwd = 2.5, xlab = 'ret arg', ylab = 'ret CF 1Y', main = "return CF -- Heston model")
lines(u.ret[,1], cf.ret[,,2], type = 'l', col = 'darkblue', lwd = 2.5)

# What if the `leverage effect' was present? The resulting CF is not symmetric around 0.5
parListHestonLev <- parListHeston
parListHestonLev$P[["1"]]$rho <- parListHestonLev$Q[["1"]]$rho <- -0.8

cf.ret.lev <- affineCF(u = u.ret, params.Q = parListHestonLev$Q, params.P = NULL, t.vec = t.vec, v.0 = v.0, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.factors = 1, CGF = FALSE)

# Take the real part
cf.ret.lev <- Re(cf.ret.lev)

lines(u.ret[,1], cf.ret.lev[,,1], type = 'l', ylim = range(cf.ret.lev), col = "red", lwd = 2.5, lty = 2)
lines(u.ret[,1], cf.ret.lev[,,2], type = 'l', col = 'darkgreen', lwd = 2.5, lty = 2)

# ---- CALCULATE MOMENTS OF STOCK RETURN ----

# Knowing the MGF allows for easy calculation of moments. This package allows for semi-closed form evaluation of the derivatives of the CF with respect to its first argument, which translates to semi-analytic calculation of the first three moments of the log-stock return.

u.mom <- matrix(c(0,0), nrow =1, ncol = 2)

t.vec.mom <- seq(1e-3,1,by=0.01)

ret.mom <- affineCFderivs(u = u.mom, params.Q = parListHeston$Q, params.P = NULL, t.vec = t.vec.mom, v.0 = v.0, jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 1)

plot(t.vec.mom, ret.mom$cf.d1[1,,1], ylim = range(ret.mom$cf.d1), type = 'l', col = 'darkorange', lwd = 2.5, xlab = "time", ylab = "E[log R]", main = 'Conditional expected log returns')
lines(t.vec.mom, ret.mom$cf.d1[1,,2], col = 'darkblue', lwd = 2.5)

ret.mom.lev <- affineCFderivs(u = u.mom, params.Q = parListHestonLev$Q, params.P = NULL, t.vec = t.vec.mom, v.0 = v.0, jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 1)

# The first moment does not change with leverage

lines(t.vec.mom, ret.mom.lev$cf.d1[1,,1], type = 'l', col = "red", lwd = 2.5, lty = 2)
lines(t.vec.mom, ret.mom.lev$cf.d1[1,,2], type = 'l', col = 'yellow', lwd = 2.5, lty = 2)

# Calculate conditional variance : second moment less first moment squared
ret.var <- ret.mom$cf.d2 - ret.mom$cf.d1^2
ret.var.lev <- ret.mom.lev$cf.d2 - ret.mom.lev$cf.d1^2

plot(t.vec.mom, ret.var[1,,1], ylim = range(c(ret.var,ret.var.lev)), type = 'l', col = 'darkorange', lwd = 2.5, xlab = "time", ylab = "V[log R]", main = 'Conditional variance of log return')
lines(t.vec.mom, ret.var[1,,2], col = 'darkblue', lwd = 2.5)
lines(t.vec.mom, ret.var.lev[1,,1], type = 'l', col = "red", lwd = 2.5, lty = 2)
lines(t.vec.mom, ret.var.lev[1,,2], type = 'l', col = 'green', lwd = 2.5, lty = 2)

# Calculate conditional skewness coefficients https://en.wikipedia.org/wiki/Skewness
ret.skew <- ret.mom$cf.d3 - 3*ret.mom$cf.d1*ret.var - ret.mom$cf.d1^3
ret.skew <- ret.skew / ret.var^1.5

ret.skew.lev <- ret.mom.lev$cf.d3 - 3*ret.mom.lev$cf.d1*ret.var.lev - ret.mom.lev$cf.d1^3
ret.skew.lev <- ret.skew.lev / ret.var.lev^1.5

plot(t.vec.mom, ret.skew[1,,1], ylim = range(c(ret.skew,ret.skew.lev)), type = 'l', col = 'darkorange', lwd = 2.5, xlab = "time", ylab = "E[(logR - E[logR])^3]/V[logR]^{1.5}", main = 'Conditional skewness of log return')
lines(t.vec.mom, ret.skew[1,,2], col = 'darkblue', lwd = 2.5)
lines(t.vec.mom, ret.skew.lev[1,,1], type = 'l', col = "red", lwd = 2.5, lty = 2)
lines(t.vec.mom, ret.skew.lev[1,,2], type = 'l', col = 'green', lwd = 2.5, lty = 2)

# ---- PLOT MGF OF VOLATILITY ----

u.vol <- cbind(0,seq(-2,2,by=0.01)) # The first argument is for the stock return, pad it with zeros.

t.vol <- 1 # plotting horizon is 1 Year

v.0 <- matrix(c(1,3),nrow = 2,ncol = 1) # We will plot the MGF at factor values 1 (the long-run mean) and 3 (a very high value)

# We will plot the risk-neutral MGF, thus params.P = NULL.
cf.vol <- affineCF(u = u.vol, params.Q = parListHeston$Q, params.P = NULL, t.vec = t.vol, v.0 = v.0, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.factors = 1, CGF = FALSE)

# Take the real part
cf.vol <- Re(cf.vol)

plot(u.vol[,2], cf.vol[,,1], type = 'l', ylim = range(cf.vol), col = "darkorange", lwd = 2.5, xlab = 'vol arg', ylab = 'vol CF 1Y', main = "variance factor CF -- Heston model")
lines(u.vol[,2], cf.vol[,,2], type = 'l', col = 'darkblue', lwd = 2.5)

# ---- MOMENTS OF VOLATILITY ----

# to calculate the conditional mean of the variance factor, you will have to take numeric derivatives

u.vol.deriv <- rbind(c(0,0),c(0,1e-3)) # The first argument is for the stock return, pad it with zeros.

# For a horizon from now to 1Y
t.vec.mom <- seq(1e-3,1,by=0.01)

# When calculating numerical derivatives, you can try pushing the ODE solver to higher precision
cf.vol.deriv <- affineCF(u = u.vol.deriv, params.Q = parListHeston$Q, params.P = NULL, t.vec = t.vec.mom, v.0 = v.0, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.factors = 1, CGF = FALSE, atol = 1e-30, rtol = 1e-14)

vol.mom <- 1e3 * (cf.vol.deriv[2,,] - cf.vol.deriv[1,,])
vol.mom <- Re(vol.mom)

plot(t.vec.mom, vol.mom[,1], ylim = range(vol.mom), type = "l", col = "darkorange", lwd = 2.5, xlab = 'time', ylab = 'E[V]', main = "Conditional expectation of vol -- Heston")
lines(t.vec.mom, vol.mom[,2], col = "darkblue", lwd = 2.5)

# ---- SIMULATE PATHS ----

# Simulate a path from the model at 100 observations per day

sim.paths <- affineSimulate(paramsList = parListHeston, N.factors = 1, t.days = 252, t.freq = 1/100, freq.subdiv = 1, rng.seed = 1, jumpGeneratorPtr = getPointerToGenerator('expNormJumpTransform'), jumpTransformPtr = getPointerToJumpTransform('expNormJumpTransform')$TF, mod.type = 'standard', nrepl = 1)

# Plot
layout(t(c(1:2)))
plot(sim.paths$S.array[,1], sim.paths$S.array[,2], type = 'l', xlab = 'day', ylab = 'Stock', main = 'Heston model \n stock price simulation')
plot(sim.paths$V.array[,1], sim.paths$V.array[,2], type = 'l', xlab = 'day', ylab = 'Vol', main = 'Heston model \n volatility simulation')
