require("rriskDistributions")
# set the design prior for binary data
# theta: the hypothesized response rate
# n: the ess of theta
# half_width of the 95% CI
design.prior.bin <- function(theta, pt, lci, uci){
  fit <- get.beta.par(p = pt, q = c(lci, theta, uci),
                      show.output = FALSE, plot = FALSE)
  a <- fit[[1]]
  b <- fit[[2]]
  c(a, b)
}

fwer1 <- function(post_prob, fwer_target = 0.10, epsilon = 0.001, gam = NULL){
  if(is.null(gam)){
    FWER <- 0
    gam <- 1
    while(FWER<=fwer_target){
      gam <- gam - epsilon
      FWER <- mean(apply(post_prob>gam, 1, any))
    }
    gam <- gam + epsilon
    FWER <- mean(apply(post_prob>gam, 1, any))
  }else{
    FWER <- mean(apply(post_prob>gam, 1, any))
  }
  return(list("gam"=gam, "FWER"=FWER))
}

# p0: hypothesized response rates for SOCs
# p1: hypothesized response rates for experimental treatments under h1
gen.scenario <- function(p0, p1){
  # assume the treatment has half the hypothesized effect
  p1_low <- p0+(p1-p0)/2
  # assume the treatment is better than the hypothesized effect
  p1_high <- p1+(p1-p0)/2
  # the number of indications
  B <- length(p0)
  
  pmat1 <- pmat2 <- pmat3 <- matrix(NA, B+1, B)
  
  for(b in 1:(B+1)){
    pmat1[b, ] <- pmat2[b, ] <- pmat3[b, ] <- p0
  }
  
  for(b in 1:B){
    # create response rates same as p1
    pmat1[-1, ][b, 1:b] <- p1[1:b]
    # create response rates lower than p1
    pmat2[-1, ][b, 1:b] <- p1_low[1:b]
    # create response rates higher than p1
    pmat3[-1, ][b, 1:b] <- p1_high[1:b]
  }
  
  rbind(pmat1, pmat2[-1, ], pmat3[-1, ])
}

##x~beta(a,b), y~beta(c,d)
##P(X>y+delta)
beta.ineq <- function(a, b, c, d, delta)
{ 
  if (a <=0 | b <= 0 | c <=0 | d <= 0) 
    stop("paramters has to be positive")
  if (a <0.01 | b < 0.01 | c <0.01 | d < 0.01) 
    stop("paramters are to close to 0")
  if (delta>1) 
    stop("delta>1!")
  if (delta<0) 
    stop("delta<0!")
  
  integrand <- function(x) { dbeta(x, a, b)*pbeta(x-delta, c, d) }
  integrate(integrand, delta, 1, rel.tol=1e-4)$value
}

## simulate data for a randomized basket trial
sim.dat <- function(n_soc, n_trt, p_soc, p_trt){
  ## sample size of each indication
  ind_size <- n_soc + n_trt
  ## number of indications
  B <- length(p_soc)
  ## response rate (trt*indication)
  pmat <- rbind(p_soc, p_trt)
  ## sample size (trt*indication)
  rownames(pmat) <- c("soc", "trt")
  colnames(pmat) <- paste("IND-", LETTERS[1:B], sep = "")
  
  rawdat <- list()
  ## aggregated data
  y_soc <- y_trt <- NULL
  for(b in 1:B){
    ## randomize within each indication
    trt <- complete_ra(N = ind_size[b],  m = n_trt[b])
    ## response probability of each patient in indication b
    pvec <- pmat[1+trt, b]
    outcome <- rbinom(n = ind_size[b], size = 1, prob = pvec)
    rawdat[[b]] <- data.frame("trt"=trt, "outcome"=outcome)
    y_trt[b] <- sum(outcome[trt==1])
    y_soc[b] <- sum(outcome[trt==0])
  }
  
  summary_data <- data.frame(n_soc = n_soc, n_trt = n_trt, 
                             y_soc = y_soc, y_trt = y_trt)
  
  list("raw_data" = rawdat, "summary_data" = summary_data)
}

#eva.trt <- function(fit, dat, a_ctr, b_ctr, B, prior_n, 
#                    method = bma_method){
#  pp <- NULL
#  for(b in 1:B){
#    beta_par <- cbind(fit$w[, b], fit$a[, b], fit$b[, b])

#    if(method == "rob"){
#      trt_fit <- mixbeta(beta_par[1, ], beta_par[2, ], beta_par[3, ])
#    }else{
#     trt_fit <- mixbeta(beta_par[1, ], beta_par[2, ])
#   }

#    w <- c(prior_n/(a_ctr+b_ctr))[b]
#    mem_prior <- c(w, 1-w)
#    mem_fit0 <- set.mem.prior(num_study = 2, delta = 0)

#    mem_fit <- update.part.bin(x = c(a_ctr[b], dat$x[1, b]), 
#                           n = c(a_ctr[b]+b_ctr[b], dat$n[1, b]), 
#                           prior_part = mem_prior, 
#                           part = mem_fit0$part)

#ctr_prior <- mixbeta(c(0.5, 1, 1), c(0.5, a_ctr[b], b_ctr[b]))
#ctr_fit <- postmix(ctr_prior, n = dat$n[1, b], r = dat$x[1, b])
#    ctr_fit <- mixbeta(c(1, 1+mem_fit$post_sim*a_ctr[b]+dat$x[1, b], 
#                   1+mem_fit$post_sim*b_ctr[b]+dat$y[1, b]))
#    pp[b] <- pmixdiff(trt_fit, ctr_fit, 0, FALSE)
#  }
# pp
#}

no.borrow <- function(dat, ni_mar){
  pp <- rep(NA, nrow(dat))
  for(b in 1:nrow(dat)){
    yt <- dat$y_trt[b]
    yc <- dat$y_soc[b]
    nt <- dat$n_trt[b]
    nc <- dat$n_soc[b]
    
    pp[b] <- beta.ineq(a = yt + 1, b = nt-yt + 1, 
                       c = yc + 1, d = nc-yc + 1, 
                       delta = ni_mar )
  }
  pp
}

## post_prob: a vector of posterior probs
solve.gam <- function(post_prob, err1_target, epsilon = 0.001){
  err1 <- 0
  gam <- 1
  while(err1 <= err1_target){
    gam <- gam - epsilon
    err1 <- mean(post_prob>gam)
  }
  gam <- gam + epsilon
  err1 <- mean(post_prob>gam)
  return(list("gam" = gam, "err1" = err1))
}

cal.pow <- function(pp_mat, gam){
  pow <- NULL
  for(i in 1:ncol(pp_mat)){
    pow[i] <- mean(pp_mat[, i]>gam[i])
  }
  pow
}

## pp_mat: matrix of posterior probs (#simulations * #indication)
## gam: posterior thereshold of different indicatons
## true_neg: indicator of true negative (TRUE vs. False) of different indications
cal.fwer <- function(pp_mat, gam, true_neg){
  new_mat <- pp_mat[, true_neg]
  new_gam <- gam[true_neg]
  if(sum(true_neg)==0){
    FWER <- 0
  }
  if(sum(true_neg)==1){
    FWER <- mean(new_mat > new_gam)
  }
  if(sum(true_neg)>1){
    FWER <- mean(apply(new_mat, 1, function(x){any(x>new_gam)}))
  }
  FWER
}