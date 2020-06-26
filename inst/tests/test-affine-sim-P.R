# test-affine simulator

library(affineModelR)

load("data/heston.params.RData")

parListHestonJmp$P$jmp$lprop <- 1
parListHestonJmp$Q$jmp$lprop <- 1
parListHestonJmp$P$jmp$muSc <- 0.0
parListHestonJmp$Q$jmp$muSc <- 0.0
parListHestonJmp$Q$`1`$eta <- 1
parListHestonJmp$Q$`1`$phi <- parListHestonJmp$P$`1`$phi
parListHestonJmp$Q$`1`$kpp <- 1.5
parListHestonJmp$P$`1`$lmb <- parListHestonJmp$Q$`1`$lmb <- 0.6
parListHestonJmp$P$`1`$rho <- parListHestonJmp$Q$`1`$rho <- -0.5

parListHestonJmp$P$`1`$erp <- 0.1
parListHestonJmp$P$`1`$erp0 <- 0.3

cf.1 <- affineCFderivs(u = matrix(0,ncol=2,nrow=1), params.Q = parListHestonJmp$Q, params.P = parListHestonJmp$P, t.vec = 21/252, v.0 = matrix(1), jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 1, mod.type = 'standard', atol = 1e-28, rtol = 1e-12)
cf.1 <- unlist(cf.1)

cf.1q <- affineCFderivs(u = matrix(0,ncol=2,nrow=1), params.Q = parListHestonJmp$Q, params.P = parListHestonJmp$Z, t.vec = 21/252, v.0 = matrix(1), jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 1, mod.type = 'standard', atol = 1e-28, rtol = 1e-12)
cf.1q <- unlist(cf.1q)

sim.1 <- affineSimulate(paramsList = parListHestonJmp, N.factors = 1, t.days = 21, t.freq = 1/78, freq.subdiv = 10, init.vals = list(S.array=0,V.array=matrix(1), day.offset=0), rng.seed = as.integer(Sys.time()), jumpGeneratorPtr = getPointerToGenerator('expNormJumpTransform'), jumpTransformPtr = getPointerToJumpTransform('expNormJumpTransform')$TF, mod.type = 'standard', nrepl = 2.5e3)

sim.1 <- lapply(sim.1$sim.arrays, function(x) x$S.array)
sim.1 <- do.call(cbind, sim.1)
sim.1 <- tail(sim.1,1)

m.test.1 <- t.test(log(sim.1), mu=cf.1[2])

v.test.1 <- pchisq(q = (length(sim.1))*((1/length(sim.1) * sum((log(sim.1) - cf.1[2])^2))/(cf.1[3] - cf.1[2]^2)), df = length(sim.1))

#### ---- Bates ----
data("affine-parameters-bates2006")

# par.bates.svj1$P$`1`$erp <- par.bates.svj1$Q$`1`$erp <- 0.0
par.bates.svj1$P$`1`$erp <- 0
par.bates.svj1$P$`1`$erp0 <- 0.15
# par.bates.svj1$P$`1`$erp0 <- par.bates.svj1$Q$`1`$erp0 <- 0

cf.b <- affineCFderivs(u = matrix(0,ncol=2,nrow=1), params.Q = par.bates.svj1$Q, params.P = par.bates.svj1$P, t.vec = 21/252, v.0 = matrix(1e-2), jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 1, mod.type = 'standard', atol = 1e-28, rtol = 1e-12)
cf.b <- unlist(cf.b)

cf.bq <- affineCFderivs(u = matrix(0,ncol=2,nrow=1), params.Q = par.bates.svj1$Q, params.P = NULL, t.vec = 21/252, v.0 = matrix(1e-2), jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 1, mod.type = 'standard', atol = 1e-28, rtol = 1e-12)
cf.bq <- unlist(cf.bq)

sim.b <- affineSimulate(paramsList = par.bates.svj1, N.factors = 1, t.days = 21, t.freq = 1/78, freq.subdiv = 10, init.vals = list(S.array=0,V.array=matrix(1e-2), day.offset=0), rng.seed = as.integer(Sys.time()), jumpGeneratorPtr = getPointerToGenerator('expNormJumpTransform'), jumpTransformPtr = getPointerToJumpTransform('expNormJumpTransform')$TF, mod.type = 'standard', nrepl = 1e4)

sim.bq <- affineSimulate(paramsList = list(P = NULL, Q = par.bates.svj1$Q), N.factors = 1, t.days = 21, t.freq = 1/78, freq.subdiv = 10, init.vals = list(S.array=0,V.array=matrix(1e-2), day.offset=0), rng.seed = as.integer(Sys.time()), jumpGeneratorPtr = getPointerToGenerator('expNormJumpTransform'), jumpTransformPtr = getPointerToJumpTransform('expNormJumpTransform')$TF, mod.type = 'standard', nrepl = 1e4)

sim.b <- lapply(sim.b$sim.arrays, function(x) x$S.array)
sim.b <- do.call(cbind, sim.b)
sim.b <- tail(sim.b,1)

m.test.b <- t.test(log(sim.b), mu=cf.b[2])

v.test.b <- pchisq(q = (length(sim.b))*((1/length(sim.b) * sum((log(sim.b) - cf.b[2])^2))/(cf.b[3] - cf.b[2]^2)), df = length(sim.b))

sim.bq <- lapply(sim.bq$sim.arrays, function(x) x$S.array)
sim.bq <- do.call(cbind, sim.bq)
sim.bq <- tail(sim.bq,1)

m.test.bq <- t.test(log(sim.bq), mu=cf.bq[2])

v.test.bq <- pchisq(q = (length(sim.bq))*((1/length(sim.bq) * sum((log(sim.bq) - cf.bq[2])^2))/(cf.bq[3] - cf.bq[2]^2)), df = length(sim.bq))

#### ---- bates and Heston ----
par.mix <- par.bates.svj1
par.mix$P[["2"]] <- parListHestonJmp$P$`1`
par.mix$Q[["2"]] <- parListHestonJmp$Q$`1`

par.mix$P[["1"]]$erp <- -0.15
par.mix$P[["1"]]$erp0 <- 0.1

par.mix$Q[["1"]]$erp <- 0
par.mix$Q[["1"]]$erp0 <- 0

par.mix$P[["2"]]$erp <- 0.05
par.mix$P[["2"]]$erp0 <- 0.0

par.mix$Q[["2"]]$erp <- 0
par.mix$Q[["2"]]$erp0 <- 0

par.mix$P[["jmp"]]$lprop <- c(20,30)
par.mix$Q[["jmp"]]$lprop <- c(20,30)

# par.mix$P$`1`$kpp <- 6
# par.mix$P$`1`$eta <- 0.012
# 
# par.mix$Q$`2`$kpp <- 1
# par.mix$Q$`2`$eta <- 1.2

par.mix$Q$jmp$sigmaYc <- 0.05
par.mix$Q$jmp$muYc <- -0.07

cf.m <- affineCFderivs(u = matrix(0,ncol=3,nrow=1), params.Q = par.mix$Q, params.P = par.mix$P, t.vec = 21/252, v.0 = matrix(c(1,1),ncol=2), jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 2, mod.type = 'standard', atol = 1e-28, rtol = 1e-12)
cf.m <- unlist(cf.m)

cf.mq <- affineCFderivs(u = matrix(0,ncol=3,nrow=1), params.Q = par.mix$Q, params.P = NULL, t.vec = 21/252, v.0 = matrix(c(1,1),ncol=2), jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 2, mod.type = 'standard', atol = 1e-28, rtol = 1e-12)
cf.mq <- unlist(cf.1q)

sim.m <- affineSimulate(paramsList = par.mix, N.factors = 2, t.days = 21, t.freq = 1/78, freq.subdiv = 10, init.vals = list(S.array=0,V.array=matrix(1,ncol=2), day.offset=0), rng.seed = as.integer(Sys.time()), jumpGeneratorPtr = getPointerToGenerator('expNormJumpTransform'), jumpTransformPtr = getPointerToJumpTransform('expNormJumpTransform')$TF, mod.type = 'standard', nrepl = 2.5e3)

sim.m <- lapply(sim.m$sim.arrays, function(x) x$S.array)
sim.m <- do.call(cbind, sim.m)
sim.m <- tail(sim.m,1)

m.test.m <- t.test(log(sim.m), mu=cf.m[2])

v.test.m <- pchisq(q = (length(sim.m))*((1/length(sim.m) * sum((log(sim.m) - cf.m[2])^2))/(cf.m[3] - cf.m[2]^2)), df = length(sim.m))


####

parListHestonJmp$P$`1`$erp0 <- 0.05

cf.2 <- affineCFderivs(u = matrix(0,ncol=2,nrow=1), params.Q = parListHestonJmp$Q, params.P = parListHestonJmp$P, t.vec = 21/252, v.0 = matrix(1), jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 1, mod.type = 'standard', atol = 1e-28, rtol = 1e-12)
cf.2 <- unlist(cf.2)

sim.2 <- affineSimulate(paramsList = parListHestonJmp, N.factors = 1, t.days = 21, t.freq = 1/78, freq.subdiv = 5, init.vals = list(S.array=0,V.array=matrix(1), day.offset=0), rng.seed = 12345, jumpGeneratorPtr = getPointerToGenerator('expNormJumpTransform'), jumpTransformPtr = getPointerToJumpTransform('expNormJumpTransform')$TF, mod.type = 'standard', nrepl = 5e3)

sim.2 <- lapply(sim.2$sim.arrays, function(x) x$S.array)
sim.2 <- do.call(cbind, sim.2)
sim.2 <- tail(sim.2,1)

m.test.2 <- t.test(log(sim.2), mu=cf.2[2])

v.test.2 <- pchisq(q = (length(sim.2))*((1/length(sim.2) * sum((log(sim.2) - cf.2[2])^2))/(cf.2[3] - cf.2[2]^2)), df = length(sim.2))


####

parListHestonJmp$P$`1`$erp0 <- 0
parListHestonJmp$P$`1`$erp <- -2

cf.3 <- affineCFderivs(u = matrix(0,ncol=2,nrow=1), params.Q = parListHestonJmp$Q, params.P = parListHestonJmp$P, t.vec = 21/252, v.0 = matrix(1), jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 1, mod.type = 'standard', atol = 1e-28, rtol = 1e-12)
cf.3 <- unlist(cf.3)

sim.3 <- affineSimulate(paramsList = parListHestonJmp, N.factors = 1, t.days = 21, t.freq = 1/78, freq.subdiv = 5, init.vals = list(S.array=0,V.array=matrix(1), day.offset=0), rng.seed = 12345, jumpGeneratorPtr = getPointerToGenerator('expNormJumpTransform'), jumpTransformPtr = getPointerToJumpTransform('expNormJumpTransform')$TF, mod.type = 'standard', nrepl = 5e3)$S.array

sim.3 <- sim.3[nrow(sim.3),seq(2,ncol(sim.3),by=2)]

m.test.3 <- t.test(log(sim.3), mu=cf.3[2])

v.test.3 <- pchisq(q = (length(sim.3))*((1/length(sim.3) * sum((log(sim.3) - cf.3[2])^2))/(cf.3[3] - cf.3[2]^2)), df = length(sim.3))

####

parListHestonJmp$P$`1`$erp0 <- -0.2
parListHestonJmp$P$`1`$erp <- -4

cf.4 <- affineCFderivs(u = matrix(0,ncol=2,nrow=1), params.Q = parListHestonJmp$Q, params.P = parListHestonJmp$P, t.vec = 21/252, v.0 = matrix(1), jumpTransform = getPointerToJumpTransform('expNormJumpTransform'), N.factors = 1, mod.type = 'standard', atol = 1e-28, rtol = 1e-12)
cf.4 <- unlist(cf.4)

sim.4 <- affineSimulate(paramsList = parListHestonJmp, N.factors = 1, t.days = 21, t.freq = 1/78, freq.subdiv = 5, init.vals = list(S.array=0,V.array=matrix(1), day.offset=0), rng.seed = 12345, jumpGeneratorPtr = getPointerToGenerator('expNormJumpTransform'), jumpTransformPtr = getPointerToJumpTransform('expNormJumpTransform')$TF, mod.type = 'standard', nrepl = 5e3)$S.array

sim.4 <- sim.4[nrow(sim.4),seq(2,ncol(sim.4),by=2)]

m.test.4 <- t.test(log(sim.4), mu=cf.4[2])

v.test.4 <- pchisq(q = (length(sim.4))*((1/length(sim.4) * sum((log(sim.4) - cf.4[2])^2))/(cf.4[3] - cf.4[2]^2)), df = length(sim.4))
