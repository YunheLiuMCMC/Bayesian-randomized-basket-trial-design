library(coda)
library(R2WinBUGS)
library(rjags)
library(R2jags)
library(dplyr)
set.seed(20202109)
source("./exnex utility.R")

####################################################
##JAGS code for the flexible exchangeable model:
##p.exch = 1, tau =0 => Complete pooled analysis
##p.exch = 1 => exchangeable or BHM
##p.exch = 0 => Complete stratified analysis
####################################################

###############################
#Binary Endpoint
###############################
Bin_rand_EXEXWB <- function(){
  
  ## EXNEX structure for treatment effect
  
  # priors for mu and tau(s)
  prec.mu.t <- pow(Prior.mu.t[2],-2)
  mu.t ~ dnorm(Prior.mu.t[1],prec.mu.t)
  
  for (j in 1:Ntau) {
    tau.prec.t[j] <- pow(Prior.tau.t[j,2],-2)
    tau.t[j] ~ dnorm(Prior.tau.t[j,1],tau.prec.t[j]) %_% I(0, )
  }
  
  # parameter model
  for (j in 1:Ngroups){
    Z.t[j] ~ dbin(p.exch.t[j],1)
    
    # theta from exch
    re.t[j] ~ dnorm(0,1)
    theta.ex.t[j] <- mu.t+re.t[j]*tau.t[tau.index[j]]
    
    # theta from non-exch
    prec.nex.t[j] <- pow(Prior.nex.t[j,2],-2)
    theta.nex.t[j]~ dnorm(Prior.nex.t[j,1], prec.nex.t[j])
    
    # select theta from latent variable
    theta.t[j] <- Z.t[j]*theta.ex.t[j]+(1-Z.t[j])*theta.nex.t[j]
  }
  
  # likelihood
  for (j in 1:Ngroups){
    y.t[j] ~ dbin(p.t[j], n.t[j])
    logit(p.t[j]) <- theta.t[j]
  }
}
#####################
#R-wrapper function 
#####################

Model_EXEX_rand <- function(y.t = NULL, n.t = NULL,
         group       = NULL, 
         group.nm    = NULL,
         Prior.mu.t  = NULL,
         Prior.tau.t = NULL,
         p.exch.t    = NULL,
         tau.index   = NULL,
         Prior.nex.t = NULL,
         pars        = c("p.t", "tau.t", "Z.t"),
         MCMC        = c(3000, 7000, 4, 1),
         DIC         = TRUE,
         seed        = 19978
)
{
  #Data
  y.t <- y.t
  n.t <- n.t
  
  Ngroups <- length(unique(group))
  Ntau <- max(tau.index)
  
  #Priors
  ##EX parts
  
  Prior.mu.t = Prior.mu.t
  
  tau.index <- if(length(tau.index)==1) rep(1, Ngroups) else tau.index
  
  if(length(tau.index) != Ngroups){stop("tau.index must be same length as number of groups")}
  
  Prior.tau.t <- rbind(NULL,  Prior.tau.t)
  
  p.exch.t <- if(length(p.exch.t)==1) rep(p.exch.t, Ngroups) else p.exch.t
  
  if(length(p.exch.t) != Ngroups){stop("p.exch.t must be same length as number of groups")}
  
  
  #NEX part
  
  Prior.nex.t <- rbind(NULL, Prior.nex.t)
  
  
  Prior.nex.t <- if(nrow(Prior.nex.t) == 1) t(replicate(Ngroups, Prior.nex.t[1,])) else Prior.nex.t
  
  #JAGS call 
  model <- Bin_rand_EXEXWB 
  
  dat <- list("y.t" = y.t, "n.t" = n.t, 
              "Prior.mu.t" = Prior.mu.t, 
              "Prior.tau.t" = Prior.tau.t,
              "tau.index" = tau.index, "p.exch.t" = p.exch.t, 
              "Prior.nex.t" = Prior.nex.t, 
              "Ngroups" = Ngroups, "Ntau" = Ntau
  )
  
  initsfun = function(i){
    list(tau = rexp(Ntau, 1))
  }
  
  inits <- lapply(rep(1, MCMC[3]),  initsfun )
  
  fit.1 <- jags(
    data = dat,
    model = model,
    parameters.to.save = pars,
    n.burnin = MCMC[1],
    n.iter = MCMC[2],
    n.chains = MCMC[3],
    n.thin = MCMC[4],
    DIC = DIC,
    jags.seed = seed, 
    quiet = TRUE
  )
  
  #Collecting samples and summary
  
  sim.1 = sim_jagswb(pars,fit.1$BUGSoutput$sims.matrix)
  sum.1 = sum_jagswb(pars,fit.1$BUGSoutput$summary)
  
  # Convergence diagnostics
  #Rhat = fit.1$BUGSoutput$summary[,"Rhat"]
  #Rhat =Rhat[grep("theta", names(Rhat))]
  
  #if(max(Rhat) > 1.1)warning("At least one Rhat Statistics in > 1.1")
  
  return(list(fit=fit.1, sim=sim.1, sum=sum.1))
  
}