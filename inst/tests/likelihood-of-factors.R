
params <- model_parameters_bates_2000_svjdc2_corrected_sigma()

params$P$`2` <- params$Q$`2` <-  NULL

params$P$jmp$lprop <- params$Q$jmp$lprop <- params$P$jmp$lprop[1]

dt <- 5/252

stock_vol_sample <- time_series_sample(model_parameters = model_parameters_bates_2000_svjdc2_corrected_sigma(), sample_length_in_days = 25100)

vmat <- stock_vol_sample$volatility_1 %>% matrix(ncol=1)
vmat_diff <- tail(vmat,-1) - head(vmat,-1)

likelihood <- function(par){
  params_0 <- params
  params_0$P$`1`$eta <- par[1]
  params_0$P$`1`$kpp <- par[2]
  params_0$P$`1`$lmb <- par[3]
  params_0$Q$`1`$lmb <- par[3]

  mod_dyn <- modelDynamics(params.P = params_0$P, params.Q = params_0$Q, dT = dt, N.factors = 1, N.points = 2, jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, mod.type = 'standard')

  means <- t(apply(vmat, 1, meanVecFun, meanListS = mod_dyn$mean.vec))
  vars <- apply(vmat, MARGIN = 1, FUN = covMatFun, covListS = mod_dyn$cov.list, covListDim = c(2,2))
  vars <- array(vars, dim = c(2,2,nrow(vmat)))
  lik <- sapply(1:(nrow(means)-1), FUN = function(ind) dnorm(x = vmat_diff[ind], mean = means[ind,2], sd = sqrt(vars[-1,-1,ind]), log = TRUE))

  sum(-lik)
}

opt <- optim(par = c(0.01484848,3.21,0.24), fn = likelihood, lower = c(0,0,0), upper = c(0.1,16,0.5), method = "L-BFGS-B")
