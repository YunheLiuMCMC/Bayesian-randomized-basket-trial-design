rm(list=ls(all=TRUE))
library(RBesT)
library(randomizr)
library(ggplot2)
library(reshape2)
library(GGally)
library(gridExtra)
library(foreach)
#library(bhmbasket)
## approxima
#library(cubature)
#library(arm)
require(LearnBayes)
#library(rjags)
library(mvtnorm)
library(xtable)
#Calculate logit transformation
#library(car)

setwd("/Users/iad/Desktop/Study/Yale/Research/Wei Wei Project TDB /code")

source("reg_bma_2.R")
source("myfun.R")
source("mem functions.R")
source("exnex.R")

#number of simulations
nsim <- 2000
r <- 0
# Noninferiority margin
ni_mar <- 0.05
## ----------------------------------------------------------------
## current trial:
## Participants will be randomized in a 2:1 ratio to an 
## experimental treatment or standard of care (SOC) based on 
## the confirmed or suspected fungal species. Below are the 
## details of comparator and sample size for each pathogen. 
##Fungal Pathogen	SOC	N (Trt/Control)
## Fusarium	L-AmB +/- VORI	20/10
## Scedosporium	VORI +/- other AF	15/7
## Lomentospora	VORI + TERB  (+/- other AF)	15/7
## Aspergillus	L-AmB +/-CAN	30/15
## Mucorales	L-AmB +/- other AF or ISA monoRx	20/10
## ----------------------------------------------------------------
# randomization ratio of trt vs. ctr
R <- 2
ind_name <- c("Fusarium", "Scedosporium", "Lomentospora", "Aspergillus", 
              "Mucorales")
# below are the assumed true death rates we simulate from:
# day42 death rates for SOC
p0 <- c(0.60, 0.88, 0.70, 0.51, 0.67)
# day42 death rates for the treatment:
p1 <- p0 - 0.30

# B: the number of indications 
B <- length(p0)
# ind_size: the total sample size of each indication
ind_size <- c(30, 22, 22, 45, 30)
n_soc <- c(10, 7, 7, 15, 10)
n_trt <- c(20, 15, 15, 30, 20)
## ----------------------------------------------------------------
## historical control:
## ----------------------------------------------------------------
## sample size of relevant historical controls:
nh <- c(112, 89, 20, 37, 21)
## number of day42 deaths in historical controls
rh <- c(66, 78, 14, 19, 14)
## prior ess for borrowing from historical trial
prior_n <- n_soc

xtable(data.frame(ind_name, rh, nh, rh/nh, n_soc, n_trt), 
       digits = c(rep(0,4), 2, rep(0, 2)))

ph <- rh/nh
## unit information of the log odds for SOC
u <- ph*(1-ph)
## sd of the prior for log odds
sd0 <- sqrt((u*n_soc)^(-1))
m0 <- log(ph/(1-ph))

par(mfrow= c(1, B))
paxis <- seq(0, 1, 0.001)
for(b in 1:B){
  pvec <- rbeta(1000, rh[b], nh[b]-rh[b])
  qvec <- log(pvec/(1-pvec))
  d <- 0.5*dnorm(seq(-5, 5, 0.01), m0[b], sd0[b])+
    0.5*dnorm(seq(-5, 5, 0.01), 0, 100)
  #plot(density(qvec), lwd = 2, col = 1, xlim = c(-2, 4), main = "", las = 1)
  plot(seq(-5, 5, 0.01), d, type = "l", lwd = 2, col = 1, xlim = c(-2, 4), main = "", las = 1)
}
## ----------------------------------------------------------------
## enthusiastic and pessimistic priors:
## design priors represent the state of knowledge about
## SOC or treatment, they should not overlap too much
## the selection of design priors can be informed by historical trials or 
## based on expert opinions. If there are multiple historical trials, 
## map can be used. 
## previous trials and external data might not have the biomarkers needed
## ----------------------------------------------
pess_prior <- enth_prior <- matrix(NA, nrow = 2, ncol = length(p1))

for(b in 1:B){
  n11 <- ind_size[b]*(R/(R+1))*p1[b]
  n10 <- ind_size[b]*(R/(R+1))*(1-p1[b])
  n01 <- ind_size[b]*(1/(R+1))*(p0[b])
  n00 <- ind_size[b]*(1/(R+1))*(1-p0[b])
  v <- sum(c(n11, n10, n01, n00)^(-1))
  #pess_prior[, b] <- c(0, sqrt(v))
  #assume it is worse
  pess_prior[, b] <- c(logit(p0[b]+ni_mar)-logit(p0[b]), sqrt(v))
  enth_prior[, b] <- c(logit(p1[b])-logit(p0[b]), sqrt(v))
}

pdf("design priors.pdf", width = 10, height = 8)
setwd("/Users/iad/Desktop/Study/Yale/Research/Wei Wei Project TDB /Paper/Paper.Plot/Plot_NI_Test_0.1_Extra")
par(mfrow=c(2, 3), oma=c(2, 4, 0.1, 1), mar=c(2, 2, 4, 2))
for(b in 1:B){
  den0 <- dnorm(seq(-6, 6, 0.01), pess_prior[1, b], pess_prior[2, b])
  den1 <- dnorm(seq(-6, 6, 0.01), enth_prior[1, b], enth_prior[2, b])
  
  plot(seq(-6, 6, 0.01), den0, type = "l", las = 1, lwd = 2, 
       xlab = "", ylab = "",  main = ind_name[b], 
       ylim = c(0, max(den0)+0.1), cex.axis = 1.5)
  lines(seq(-6, 6, 0.01), den1, col = 2, lwd = 2)
  
  if(b==1){
    title(xlab = "", ylab = "Prior Density", cex.lab= 1.5, 
          outer = T, line = 2)
    legend("topright", lty = 1, lwd = 2, col = 1:2, 
           c("Pessimistic", "Enthusiastic"), bty = "n", cex = 1.5)
  }
  #if(b==B){
  #  title(xlab = "Log Odds Ratio", ylab = "", cex.lab= 1.5, 
  #       outer = T, line = 0.6)
  #}
}
dev.off()
## ----------------------------------------------------------------
## set up the simulation study
## ----------------------------------------------------------------
# pmat0/pmat1: matrix of response rates for the SOC/trt 
## rows: scenarios
## columns: indications
#pmat0 <- rbind(gen.scenario(p0 = p0, p1 = p0))[1:(B+1), ]
#pmat1 <- rbind(gen.scenario(p0 = p0, p1 = p1))[1:(B+1), ]
pmat0 <- rbind(gen.scenario(p0 = p0, p1 = p0))
pmat1 <- rbind(gen.scenario(p0 = p0, p1 = p1))
##non-inferiority test
pmat1[1, ] <- pmat1[1, ] + 0.05
for(i in 2:5){
  pmat1[i, i:B] <- pmat1[i, i:B] + 0.05
}
for(i in 2:5){
  pmat1[i+5, i:B] <- pmat1[i+5, i:B] + 0.05
}
for(i in 2:5){
  pmat1[i+10, i:B] <- pmat1[i+10, i:B] + 0.05
}


#S: the number of scenarios
pmat0 <- rbind(pmat0[1:6, ], c(0.6, 0.88, 0.7, 0.51, 0.67), pmat0[7:11, ])
pmat1 <- rbind(pmat1[1:6, ], c(0.6, 0.88, 0.7, 0.51, 0.67), pmat1[7:11, ])
S <- nrow(pmat0)

## ----------------------------------------------------------------
## select the method for BMA
## ----------------------------------------------------------------
bma_mods <- init.mod(B, prior_indicator = c(0, 1))

bma_prior <- gen.prior(mods = bma_mods, r = r)$prior

aggregate(bma_prior, by = list(bma_mods[, 1]), sum)
aggregate(bma_prior, by = list(rowSums(bma_mods)), sum)

mat <- cbind(gen.prior(mods = bma_mods, r = 0)$prior, 
             gen.prior(mods = bma_mods, r = 1)$prior, 
             gen.prior(mods = bma_mods, r = 2)$prior)

sim0 <- gen.prior(mods = bma_mods, r = 0)$sim
sim1 <- gen.prior(mods = bma_mods, r = 1)$sim
sim2 <- gen.prior(mods = bma_mods, r = 2)$sim

g0 <- ggcorr(data = NULL, cor_matrix = sim0, nbreaks=5, label=TRUE, 
             label_size=5, size=4, legend.size=5, label_round = 2, 
             legend.position = "none")

g1 <- ggcorr(data = NULL, cor_matrix = sim1, nbreaks=5, label=TRUE, 
             label_size=5, size=4, legend.size=10, label_round = 2, 
             legend.position = "none")

g2 <- ggcorr(data = NULL, cor_matrix = sim2, nbreaks=5, label=TRUE, 
             label_size=5, size=4, legend.size=10, label_round = 2, 
             legend.position = "none")

fig <- grid.arrange(g0, g1, g2, nrow = 1)
#ggsave(file="prior_corr.pdf", fig)

#logistic.bma(dat = tmp, b0_prior <- c(0, 5, 1), pess_prior = pess_prior, 
 #            enth_prior = enth_prior, mod_prior = bma_prior)$mlik

