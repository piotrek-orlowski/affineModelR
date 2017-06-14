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

mkt <- data.frame(t=t.vec,r=0.03,q=0,p=1)
mkt.0 <- data.frame(t=t.vec,r=0.0,q=0,p=1)

v.0 <- matrix(c(1,3),nrow = 2,ncol = 1) # We will plot the MGF at factor values 1 (the long-run mean) and 3 (a very high value)

# We will plot the risk-neutral MGF, thus params.P = NULL.
cf.ret <- affineCF(u = u.ret, params.Q = parListHeston$Q, params.P = NULL, t.vec = NULL, v.0 = v.0, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.factors = 1, CGF = FALSE, mkt = mkt)

# Take the real part
cf.ret <- Re(cf.ret)

plot(u.ret[,1], cf.ret[,,1], type = 'l', ylim = range(cf.ret), col = "darkorange", lwd = 2.5, xlab = 'ret arg', ylab = 'ret CF 1Y', main = "return CF -- Heston model")

cf.der <- affineCFderivs(u = u.ret, params.Q = parListHeston$Q, params.P = NULL, t.vec = NULL, v.0 = v.0, jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 1, CGF = FALSE, mkt = mkt)

u.0 <- cbind(0,0)

# with mkt.0
cf.d0 <- affineCFderivs(u = u.0, params.Q = parListHeston$Q, params.P = NULL, t.vec = NULL, v.0 = v.0, jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 1, CGF = FALSE, mkt = mkt.0, atol = 1e-37, rtol = 1e-15)

cf.d0n <- affineCFderivsNumerical(u = u.0, params.Q = parListHeston$Q, params.P = NULL, t.vec = NULL, v.0 = v.0, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.factors = 1, mkt = mkt.0, mod.type = 'standard', atol = 1e-37, rtol = 1e-15)

# with mkt, i.e. interest rate
cf.drn <- affineCFderivsNumerical(u = u.0, params.Q = parListHeston$Q, params.P = NULL, t.vec = NULL, v.0 = v.0, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.factors = 1, mkt = mkt, mod.type = 'standard', atol = 1e-37, rtol = 1e-15)

# first derivative with interest rates
dim(cf.d0$cf.d1)
ermqt <- exp(-mkt$t*(mkt$r-mkt$q))
rmqt <- mkt$t*(mkt$r-mkt$q)
df <- array(exp(rmqt*u.0[,1]),dim = dim(cf.d0$cf))

cf.d0a <- df * (cf.d0$cf.d1 + rmqt * cf.d0$cf)

# second derivative with interest rates
cf.d0a2 <- rmqt * cf.d0a + df *(cf.d0$cf.d2 + rmqt * cf.d0$cf.d1)

# third derivative with interest rates
cf.d0a3 <- rmqt * cf.d0a2 + rmqt * df *(cf.d0$cf.d2 + rmqt * cf.d0$cf.d1) + df*(cf.d0$cf.d3 + rmqt * cf.d0$cf.d2)
