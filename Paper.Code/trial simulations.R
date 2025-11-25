# array to store results from S scenarios
ref_pp <- bma_pp <- bma2D_pp <- nex_pp <- nex2D_pp <- array(NA, c(S, nsim, B))

pool_zval <- test_pool_zval <- array(NA, c(S, nsim, B))

# Matrix for saving the power of each basket in each scenario
ref_pow <- bma_pow <- bma2D_pow <- nex_pow <- nex2D_pow <- matrix(NA, nrow = S, ncol = B) 

# Informative prior for the control arm using exnex
infor_prior_ab <- cbind(prior_n * rh/nh, prior_n * (nh-rh)/nh)

# array to store bias from S scenarios
ref_bias <- pool_bias <- pool_test_bias <- bma_bias <- bma2D_bias <- nex_bias <- nex2D_bias <- array(NA, c(S, nsim, B))

## Hypothesis of NI test:
## H0: p1 >= p0 + 0.05 (the D42 mortality rate is higher than p0+0.05)
## H1: p1 < p0 + 0.05
for(s in 1:S){
  for(i in 1:nsim){
    #set.seed(i)
    sim_out <- sim.dat(n_soc = n_soc, n_trt = n_trt, 
                       p_soc = pmat0[s, ] + ni_mar, p_trt = pmat1[s, ])
    
    trial_dat <- sim_out$raw_data
    ## ----------------------------------------------------
    ## Reference method without any borrowing
    ## ----------------------------------------------------
    trial_summary <- sim_out$summary_data
    trial_summary$Disease <- c("Fusarium", "Scedosporium", "Lomentospora", "Aspergillus", "Mucorales")
    ref_pp[s, i, ] <- 1-no.borrow(trial_summary, ni_mar = ni_mar)
    #Calculate mortality rate est bias
    phat0 <- (sim_out$summary_data$y_soc+1)/(sim_out$summary_data$n_soc+2)
    phat1 <- (sim_out$summary_data$y_trt+1)/(sim_out$summary_data$n_trt+2)
    bias.ref <- (phat1-phat0) - (pmat1[s, ]-pmat0[s, ])
    ref_bias[s, i, ] <- bias.ref
    ## ----------------------------------------------------
    ## naive pooling
    ## ----------------------------------------------------
    #pooled <- data.frame(matrix(colSums(sim_out$summary_data), nrow = 1))
    #names(pooled) <- names(sim_out$summary_data)
    pooled <- colSums(sim_out$summary_data)
    phat0 <- pooled[3]/pooled[1]
    phat1 <- pooled[4]/pooled[2]
    se <- sqrt(phat0*(1-phat0)/pooled[1] + phat1*(1-phat1)/pooled[2])
    pool_zval[s, i, ] <- (phat1 - phat0 - ni_mar )/se
    #pool_pp[i, s] <- 1-no.borrow(pooled, ni_mar = ni_mar)
    #Calculate mortality rate est bias
    bias.pool <- (phat1-phat0) - (pmat1[s, ]-pmat0[s, ])
    pool_bias[s, i, ] <- bias.pool
    
    ## ----------------------------------------------------
    ## if pooled result is significant, then proceed with individual tests:
    ## ----------------------------------------------------
    x <- as.vector(sim_out$summary_data$n_trt - sim_out$summary_data$y_trt)
    y <- as.vector(sim_out$summary_data$y_trt)
    tab <- cbind(x, y)
    #pval <- fisher.test(tab)$p.val
    ##zvalue to pvalue
    pval <- pnorm(unique(pool_zval[s, i, ]))
    if(pval>0.1){
      test_pool_zval[s, i, ] <- pool_zval[s, i, ]
      #Calculate mortality rate est bias
      pool_test_bias[s, i, ] <- pool_bias[s, i, ]
    }
    if(pval<0.1){
      phat0 <- sim_out$summary_data$y_soc/sim_out$summary_data$n_soc
      phat1 <- sim_out$summary_data$y_trt/sim_out$summary_data$n_trt
      tot_n <- sim_out$summary_data$n_soc + sim_out$summary_data$n_trt
      tot_y <- sim_out$summary_data$y_soc + sim_out$summary_data$y_trt
      phat <- tot_y/tot_n
      se <- sqrt(phat0*(1-phat0)/sim_out$summary_data$n_soc + phat1*(1-phat1)/sim_out$summary_data$n_trt)
      
      test_pool_zval[s, i, ] <- (phat1 - phat0 - ni_mar)/se
      
      #Calculate mortality rate est bias
      phat0 <- sim_out$summary_data$y_soc/sim_out$summary_data$n_soc
      phat1 <- sim_out$summary_data$y_trt/sim_out$summary_data$n_trt
      bias.pool.test <- (phat1-phat0) - (pmat1[s, ]-pmat0[s, ])
      pool_test_bias[s, i, ] <- bias.pool.test
    }
    
    
    ## ----------------------------------------------------
    ## BMA: learning across indications
    ## ----------------------------------------------------
    b0_prior_mat <- cbind(rep(0, B), rep(5, B), 1)
    
    post_mod_wt <- update.mod.wt(dat = trial_dat, b0_prior_mat = b0_prior_mat, pess_prior = pess_prior, 
                                 enth_prior = enth_prior, mod_prior = bma_prior)$mods
    
    post_enth_wt <- colSums(post_mod_wt[, 1:B]*post_mod_wt$p)
    
    beta_ni <- log((p0+ni_mar)/(1-(p0+ni_mar)))-log(p0/(1-p0))
    
    fit <- fit.bma(dat = trial_dat, b0_prior_mat = b0_prior_mat, pess_prior = pess_prior, 
                   enth_prior = enth_prior, mod_prior = bma_prior,
                   post_enth_wt = post_enth_wt, beta_ni = beta_ni)
    bma_pp[s, i, ] <- fit$pp
    
    #Calculate mortality rate est bias
    bias.bma <- fit$bias - (pmat1[s, ]-pmat0[s, ])
    bma_bias[s, i, ] <- bias.bma 
    #summary(glm(tmp[[3]]$outcome ~ tmp[[3]]$trt, family = binomial))$coef
    ## ----------------------------------------------------
    ## two dimensional learning 
    ## ----------------------------------------------------
    b0_prior_mat <- cbind(m0, sd0, 0.5)
    
    post_mod_wt <- update.mod.wt(dat = trial_dat, b0_prior_mat = b0_prior_mat, pess_prior = pess_prior, 
                                 enth_prior = enth_prior, mod_prior = bma_prior)$mods
    
    post_enth_wt <- colSums(post_mod_wt[, 1:B]*post_mod_wt$p)
    
    fit <- fit.bma(dat = trial_dat, b0_prior_mat = b0_prior_mat, pess_prior = pess_prior, 
                   enth_prior = enth_prior, mod_prior = bma_prior,
                   post_enth_wt = post_enth_wt, beta_ni = beta_ni)
    
    bma2D_pp[s, i, ] <- fit$pp
    
    #Calculate mortality rate est bias
    bias.bma2D <- fit$bias - (pmat1[s, ]-pmat0[s, ])
    bma2D_bias[s, i, ] <- bias.bma2D 
    ## ----------------------------------------------------
    ## EXNEX
    ## ----------------------------------------------------
    trial_summary$x_soc <- trial_summary$n_soc-trial_summary$y_soc
    
    nex_anal <- Model_EXEX_rand (y.t = trial_summary$y_trt,
                                 n.t = trial_summary$n_trt,
                                 group       = 1:B, 
                                 group.nm    = trial_summary$Disease,
                                 Prior.mu.t  = c(log(4), 1),
                                 Prior.tau.t = c(0, 0.5),
                                 p.exch.t    = 0.5,
                                 tau.index   = 1,
                                 Prior.nex.t = cbind(rep(0, B), rep(2, B)),
                                 pars = c("p.t", "tau.t", "Z.t"),
                                 MCMC = c(5000, 10000, 1, 1),
                                 DIC =TRUE, seed = 19978)
    
    pt_mcmc <- nex_anal$sim$p.t 
    
    pc_mcmc <- NULL
    for(basket in 1:B){
      pc_mcmc <- cbind(pc_mcmc, rbeta(5000, 1+trial_summary$y_soc[basket], 1+trial_summary$x_soc[basket]))
    }
    
    nex_pp[s, i, ] <- colMeans(pt_mcmc < pc_mcmc + ni_mar)
    
    #Calculate mortality rate est bias
    bias.nex <- (colMeans(pt_mcmc) - colMeans(pc_mcmc)) - (pmat1[s, ]-pmat0[s, ])
    nex_bias[s, i, ] <- bias.nex
    
    ## ----------------------------------------------------
    ## EXNEX-2D
    ## ----------------------------------------------------
    pc_mcmc <- NULL
    for(basket in 1:B){
      ## beta mixture priors 
      bm_prior <- mixbeta(rob = c(0.5, 1, 1), inf = c(0.5, infor_prior_ab[basket, ]))
      bm_post <- postmix(priormix = bm_prior, n = trial_summary$n_soc[basket], 
                         r = trial_summary$y_soc[basket])
      
      pc_mcmc <- cbind(pc_mcmc, rmix(bm_post, 5000))
    }
    
    nex2D_pp[s, i, ] <- colMeans(pt_mcmc < pc_mcmc + ni_mar)
    
    #Calculate mortality rate est bias
    bias.nex2D <- (colMeans(pt_mcmc) - colMeans(pc_mcmc)) - (pmat1[s, ]-pmat0[s, ])
    nex2D_bias[s, i, ] <- bias.nex2D
    
    print(c(i, s))
  }
}

setwd("/Users/iad/Desktop")
if(ni_mar==0){
  save.image(paste("IMI Extra results r=", r, ".RData", sep = ""))
}else{
  save.image(paste("IMI Extra NI results r=", r, ".RData", sep = ""))
}

