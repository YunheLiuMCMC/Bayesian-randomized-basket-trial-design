rm(list=ls(all=TRUE))
library(xtable)
setwd("/Users/iad/Desktop/Study/Yale/Research/Wei Wei Project TDB /Paper/Paper.Data/NI.Test.0.1.Extra")
load("IMI Extra NI results r=0.RData")
#######################Tabulate fwer && Visualized power########################
setwd("/Users/iad/Desktop/Study/Yale/Research/Wei Wei Project TDB /Paper/Paper.Plot/Plot_NI_Test_0.1_Extra")
ref_pow <- matrix(NA, nrow = S, ncol = B)
bma_pow <- bma2D_pow <- matrix(NA, nrow = S, ncol = B)
nex_pow <- nex2D_pow <- matrix(NA, nrow = S, ncol = B)
test_pool_pow <- pool_pow <- matrix(NA, nrow = S, ncol = B)

fwer_ref <- fwer_bma <- fwer_bma2D <- NULL
fwer_nex <- fwer_nex2D <- fwer_pool <- fwer_test_pool <- NULL
  
fwer_target <- 0.1

gam_ref <- rep(fwer1(post_prob = ref_pp[1, ,], fwer_target = fwer_target, epsilon = 0.0001)$gam, B)
gam_bma <- rep(fwer1(post_prob = bma_pp[1, ,], fwer_target = fwer_target, epsilon = 0.0001)$gam, B)
gam_bma2D <- rep(fwer1(post_prob = bma2D_pp[1, ,], fwer_target = fwer_target, epsilon = 0.0001)$gam, B)
gam_nex <- rep(fwer1(post_prob =   nex_pp[1, ,],   fwer_target = fwer_target, epsilon = 0.0001)$gam, B)
gam_nex2D <- rep(fwer1(post_prob = nex2D_pp[1, ,], fwer_target = fwer_target, epsilon = 0.0001)$gam, B)


for(s in 1:S){
  ref_pow[s, ] <-   cal.pow(pp_mat = ref_pp[s, ,],   gam = gam_ref)
  bma_pow[s, ] <-   cal.pow(pp_mat = bma_pp[s, ,],   gam = gam_bma)
  bma2D_pow[s, ] <- cal.pow(pp_mat = bma2D_pp[s, ,], gam = gam_bma2D)
  nex_pow[s, ] <-   cal.pow(pp_mat = nex_pp[s, ,],   gam = gam_nex)
  nex2D_pow[s, ] <- cal.pow(pp_mat = nex2D_pp[s, ,], gam = gam_nex2D)
  
  pool_pow[s, ] <- colMeans(pool_zval[s, , ] < qnorm(fwer_target))
  
  pool <- apply(test_pool_zval[s, , ] - pool_zval[s, , ], 1, function(x){all(x==0)})
  tmp1 <- test_pool_zval[s, pool==TRUE, ] < qnorm(fwer_target)
  tmp2 <- test_pool_zval[s, pool==FALSE, ] < qnorm(fwer_target/B)
  test_pool_pow[s, ] <- colMeans(rbind(tmp1, tmp2))
  
  #NI
  idx <- pmat1[s, ]==pmat0[s, ] + ni_mar
  fwer_ref[s] <-   cal.fwer(pp_mat = ref_pp[s, , ],   gam = gam_ref,   true_neg = idx)
  fwer_bma[s] <-   cal.fwer(pp_mat = bma_pp[s, , ],   gam = gam_bma,   true_neg = idx)
  fwer_bma2D[s] <- cal.fwer(pp_mat = bma2D_pp[s, , ], gam = gam_bma2D, true_neg = idx)
  fwer_nex[s] <-   cal.fwer(pp_mat = nex_pp[s, , ],   gam = gam_nex,   true_neg = idx)
  fwer_nex2D[s] <- cal.fwer(pp_mat = nex2D_pp[s, , ], gam = gam_nex2D, true_neg = idx)
  
  if(sum(idx)==1){
    fwer_test_pool[s] <- mean(rbind(tmp1, tmp2)[, idx])
    fwer_pool[s] <- mean(pool_zval[s, , idx] < qnorm(fwer_target))
  }else{
    fwer_test_pool[s] <- mean(apply(rbind(tmp1, tmp2)[, idx], 1, any))
    
    fwer_pool[s] <- mean(apply(pool_zval[s, , idx] < qnorm(fwer_target), 1, any))
  }
  
  print(mean(pool))
}

ref_tab <- cbind(ref_pow, fwer_ref)
bma_tab <- cbind(bma_pow, fwer_bma)
bma2D_tab <- cbind(bma2D_pow, fwer_bma2D)
nex_tab <- cbind(nex_pow, fwer_nex)
nex2D_tab <- cbind(nex2D_pow, fwer_nex2D)
test_pool_tab <- cbind(test_pool_pow, fwer_test_pool)
pool_tab <- cbind(pool_pow, fwer_pool)

#ref_tab
#bma_tab
#bma2D_tab
#nex_tab
#nex2D_tab
#xtable(cbind(pmat, result$fw[, 1:4]), digits = c(0, rep(3, 8))) 
#pwr <- c(as.vector(t(ref_pow)), 
#         as.vector(t(bma_pow)), as.vector(t(bma2D_pow)), 
 #        as.vector(t(nex_pow)), as.vector(t(nex2D_pow)),
 #        as.vector(t(pool_test_pow)))
pwr <- c(as.vector(t(ref_pow)),
         as.vector(t(bma_pow)), as.vector(t(bma2D_pow)), 
         as.vector(t(nex_pow)), as.vector(t(nex2D_pow)))
met <- rep(c("REF", "BMA", "BMA-2D", "EXNEX", "EXNEX-2D"), each = S*B)
ind <- rep(rep(1:B, S), 5)
sce <- rep(rep(1:S, each = B), 5)

bar_dat <- data.frame(pwr, met, ind, sce)
bar_dat$met <- factor(bar_dat$met, levels = c("REF", "BMA", "EXNEX", "BMA-2D", "EXNEX-2D"))
power_name_1 <- paste0("NI_power_1_","r_",r,".pdf")
pdf(power_name_1, width = 12, height = 8)
# power under mixed alt
colors <- c("#8dd3c7", "#bebada", "#fb8072", "#80b1d3", "#fdb462")
Indication_Name <- c("F", "S", "L", "A", "M")
par(mfrow=c(2, 3), oma=c(1,5,0.1,2), mar=c(4,2,2,2))
for(s in 2:6){
  barplot(pwr ~ met + ind, data = bar_dat[bar_dat$sce==s, ], las = 1, 
          ylab = "", xlab = "", xaxt = "n", ylim = c(0, 1),
          beside=TRUE,  cex.axis = 1.25, 
          col = colors,
          main = paste("S", s-1, sep=""), cex.lab = 1.5)
  abline(h = 0.1, lty = 3, col = 2)
  axis(1, seq(3.5, 30, 6)[1:B], paste("IND-", Indication_Name, sep=""))
  #mtext(paste(pmat0[s, ],"(", pmat1[s, ], ")", sep=""), side=1, cex=0.7, line=2.5, at=seq(3, 25, 5)[1:B], font=2)
  if(s==2){
    legend("topright", pch = 15, cex = 2, bty = "",
           col = colors,
           c("REF", "BMA", "EXNEX", "BMA-2D", "EXNEX-2D"))
  }
}
dev.off()

power_name_2 <- paste0("NI_power_2_","r_",r,".pdf")
pdf(power_name_2, width = 12, height = 8)
# power under mixed alt
par(mfrow=c(2, 3), oma=c(1,5,0.1,2), mar=c(4,2,2,2))
for(s in 7:12){
  barplot(pwr ~ met + ind, data = bar_dat[bar_dat$sce==s, ], las = 1, 
          ylab = "", xlab = "", xaxt = "n", ylim = c(0, 1),
          beside=TRUE,  cex.axis = 1.25, 
          col = colors,
          main = paste("S", s-1, sep=""), cex.lab = 1.5)
  abline(h = 0.1, lty = 3, col = 2)
  axis(1, seq(3.5, 30, 6)[1:B], paste("IND-", Indication_Name, sep=""))
  #mtext(paste(pmat0[s, ],"(", pmat1[s, ], ")", sep=""), side=1, cex=0.7, line=2.5, at=seq(3, 25, 5)[1:B], font=2)
  if(s==7){
    legend("topright", pch = 15, cex = 2, bty = "",
           col = colors,
           c("REF", "BMA", "EXNEX", "BMA-2D", "EXNEX-2D"))
  }
}
dev.off()

#result$fw
fwer_tab <- cbind(fwer_ref, fwer_bma, fwer_bma2D, fwer_nex, fwer_nex2D, fwer_test_pool, fwer_pool)
fwer_tab
##############Summary BMA and BMA-2D Sensitivity Analysis Results###############
##Store power
ref_pow <- matrix(NA, nrow = S, ncol = B)
bma_pow <- bma2D_pow <- matrix(NA, nrow = S, ncol = B)
##Store fwer
fwer_ref <- fwer_bma <- fwer_bma2D <- NULL
fwer_nex <- fwer_nex2D <- fwer_pool <- fwer_test_pool <- NULL
##Define fwer
fwer_target <- 0.1
################################Read r=0.Rdata##################################
setwd("/Users/iad/Desktop/Study/Yale/Research/Wei Wei Project TDB /Paper/Paper.Data/NI.Test.0.1.Extra")
load("IMI Extra NI results r=0.RData")
##Calculate the tuning parameter to control fwer
gam_ref <- rep(fwer1(post_prob = ref_pp[1, ,], fwer_target = fwer_target, epsilon = 0.0001)$gam, B)
gam_bma <- rep(fwer1(post_prob = bma_pp[1, ,], fwer_target = fwer_target, epsilon = 0.0001)$gam, B)
gam_bma2D <- rep(fwer1(post_prob = bma2D_pp[1, ,], fwer_target = fwer_target, epsilon = 0.0001)$gam, B)

for(s in 1:S){
  ref_pow[s, ] <-   cal.pow(pp_mat = ref_pp[s, ,],   gam = gam_ref)
  bma_pow[s, ] <-   cal.pow(pp_mat = bma_pp[s, ,],   gam = gam_bma)
  bma2D_pow[s, ] <- cal.pow(pp_mat = bma2D_pp[s, ,], gam = gam_bma2D)
}

AVP_0 <- matrix(NA, nrow = S, ncol = 3)
for (s in 1:S){
  if (s == 1){
    AVP_0[s, 1] <- mean(ref_pow[s, ])
    AVP_0[s, 2] <- mean(bma_pow[s, ])
    AVP_0[s, 3] <- mean(bma2D_pow[s, ])
  }
  if (s >1 & s <= 6){
    AVP_0[s, 1] <- mean(ref_pow[s, 1:(s-1)])
    AVP_0[s, 2] <- mean(bma_pow[s, 1:(s-1)])
    AVP_0[s, 3] <- mean(bma2D_pow[s, 1:(s-1)])
  }
  if (s == 7){
    AVP_0[s, 1] <- mean(ref_pow[s, ])
    AVP_0[s, 2] <- mean(bma_pow[s, ])
    AVP_0[s, 3] <- mean(bma2D_pow[s, ])
  }
  if (s > 7){
    AVP_0[s, 1] <- mean(ref_pow[s, 1:(s-7)])
    AVP_0[s, 2] <- mean(bma_pow[s, 1:(s-7)])
    AVP_0[s, 3] <- mean(bma2D_pow[s, 1:(s-7)])
  }
}
################################Read r=-2.Rdata##################################
setwd("/Users/iad/Desktop/Study/Yale/Research/Wei Wei Project TDB /Paper/Paper.Data/NI.Test.0.1.Extra")
load("IMI Extra NI results r=-2.RData")
##Calculate the tuning parameter to control fwer
gam_ref <- rep(fwer1(post_prob = ref_pp[1, ,], fwer_target = fwer_target, epsilon = 0.0001)$gam, B)
gam_bma <- rep(fwer1(post_prob = bma_pp[1, ,], fwer_target = fwer_target, epsilon = 0.0001)$gam, B)
gam_bma2D <- rep(fwer1(post_prob = bma2D_pp[1, ,], fwer_target = fwer_target, epsilon = 0.0001)$gam, B)

for(s in 1:S){
  ref_pow[s, ] <-   cal.pow(pp_mat = ref_pp[s, ,],   gam = gam_ref)
  bma_pow[s, ] <-   cal.pow(pp_mat = bma_pp[s, ,],   gam = gam_bma)
  bma2D_pow[s, ] <- cal.pow(pp_mat = bma2D_pp[s, ,], gam = gam_bma2D)
}

AVP_m_2 <- matrix(NA, nrow = S, ncol = 3)
for (s in 1:S){
  if (s == 1){
    AVP_m_2[s, 1] <- mean(ref_pow[s, ])
    AVP_m_2[s, 2] <- mean(bma_pow[s, ])
    AVP_m_2[s, 3] <- mean(bma2D_pow[s, ])
  }
  if (s >1 & s <= 6){
    AVP_m_2[s, 1] <- mean(ref_pow[s, 1:(s-1)])
    AVP_m_2[s, 2] <- mean(bma_pow[s, 1:(s-1)])
    AVP_m_2[s, 3] <- mean(bma2D_pow[s, 1:(s-1)])
  }
  if (s == 7){
    AVP_m_2[s, 1] <- mean(ref_pow[s, ])
    AVP_m_2[s, 2] <- mean(bma_pow[s, ])
    AVP_m_2[s, 3] <- mean(bma2D_pow[s, ])
  }
  if (s > 7){
    AVP_m_2[s, 1] <- mean(ref_pow[s, 1:(s-7)])
    AVP_m_2[s, 2] <- mean(bma_pow[s, 1:(s-7)])
    AVP_m_2[s, 3] <- mean(bma2D_pow[s, 1:(s-7)])
  }
}
################################Read r=2.Rdata##################################
setwd("/Users/iad/Desktop/Study/Yale/Research/Wei Wei Project TDB /Paper/Paper.Data/NI.Test.0.1.Extra")
load("IMI Extra NI results r=2.RData")
##Calculate the tuning parameter to control fwer
gam_ref <- rep(fwer1(post_prob = ref_pp[1, ,], fwer_target = fwer_target, epsilon = 0.0001)$gam, B)
gam_bma <- rep(fwer1(post_prob = bma_pp[1, ,], fwer_target = fwer_target, epsilon = 0.0001)$gam, B)
gam_bma2D <- rep(fwer1(post_prob = bma2D_pp[1, ,], fwer_target = fwer_target, epsilon = 0.0001)$gam, B)

for(s in 1:S){
  ref_pow[s, ] <-   cal.pow(pp_mat = ref_pp[s, ,],   gam = gam_ref)
  bma_pow[s, ] <-   cal.pow(pp_mat = bma_pp[s, ,],   gam = gam_bma)
  bma2D_pow[s, ] <- cal.pow(pp_mat = bma2D_pp[s, ,], gam = gam_bma2D)
}

AVP_2 <- matrix(NA, nrow = S, ncol = 3)
for (s in 1:S){
  if (s == 1){
    AVP_2[s, 1] <- mean(ref_pow[s, ])
    AVP_2[s, 2] <- mean(bma_pow[s, ])
    AVP_2[s, 3] <- mean(bma2D_pow[s, ])
  }
  if (s >1 & s <= 6){
    AVP_2[s, 1] <- mean(ref_pow[s, 1:(s-1)])
    AVP_2[s, 2] <- mean(bma_pow[s, 1:(s-1)])
    AVP_2[s, 3] <- mean(bma2D_pow[s, 1:(s-1)])
  }
  if (s == 7){
    AVP_2[s, 1] <- mean(ref_pow[s, ])
    AVP_2[s, 2] <- mean(bma_pow[s, ])
    AVP_2[s, 3] <- mean(bma2D_pow[s, ])
  }
  if (s > 7){
    AVP_2[s, 1] <- mean(ref_pow[s, 1:(s-7)])
    AVP_2[s, 2] <- mean(bma_pow[s, 1:(s-7)])
    AVP_2[s, 3] <- mean(bma2D_pow[s, 1:(s-7)])
  }
}
################################################################################
AVP_0
AVP_m_2
AVP_2
#######Visualize different methods' mean absolute bias over trial level#########
rm(list=ls(all=TRUE))
setwd("/Users/iad/Desktop/Study/Yale/Research/Wei Wei Project TDB /Paper/Paper.Data/NI.Test.0.1.Extra")
load("IMI Extra NI results r=0.RData")
setwd("/Users/iad/Desktop/Study/Yale/Research/Wei Wei Project TDB /Paper/Paper.Plot/Plot_NI_Test_0.1_Extra")

ref_bias_ave <- matrix(NA, nrow = S, ncol = nsim)
pool_bias_ave <- matrix(NA, nrow = S, ncol = nsim)
pool_test_bias_ave <- matrix(NA, nrow = S, ncol = nsim)
bma_bias_ave <- matrix(NA, nrow = S, ncol = nsim)
nex_bias_ave <- matrix(NA, nrow = S, ncol = nsim)
bma2D_bias_ave <- matrix(NA, nrow = S, ncol = nsim)
nex2D_bias_ave <- matrix(NA, nrow = S, ncol = nsim)

for(s in 1:S){
  ref_bias_ave[s, ] <- rowMeans(abs(ref_bias[s, ,]))
  pool_bias_ave[s, ] <- rowMeans(abs(pool_bias[s, ,]))
  pool_test_bias_ave[s, ] <- rowMeans(abs(pool_test_bias[s, ,]))
  bma_bias_ave[s, ] <- rowMeans(abs(bma_bias[s, ,]))
  nex_bias_ave[s, ] <- rowMeans(abs(nex_bias[s, ,]))
  bma2D_bias_ave[s, ] <- rowMeans(abs(bma2D_bias[s, ,]))
  nex2D_bias_ave[s, ] <- rowMeans(abs(nex2D_bias[s, ,]))
}

bias <- c(as.vector(t(ref_bias_ave)), as.vector(t(pool_bias_ave)), as.vector(t(pool_test_bias_ave)),
          as.vector(t(bma_bias_ave)), as.vector(t(nex_bias_ave)),
          as.vector(t(bma2D_bias_ave)), as.vector(t(nex2D_bias_ave)))
met <- rep(c("REF", "Pool", "Test&Pool",  "BMA", "EXNEX", "BMA-2D", "EXNEX-2D"), each = S*nsim)
sce <- rep(rep(1:S, each = nsim), 7)

box_dat <- data.frame(bias, met, sce)
box_dat$met <- factor(box_dat$met, levels = c("REF", "Pool", "Test&Pool",  "BMA", "EXNEX", "BMA-2D", "EXNEX-2D"))

##
bias_name_1 <- paste0("NI_bias_1_","r_",r,".pdf")
pdf(bias_name_1, width = 12, height = 8)
par(mfrow=c(2, 3), oma=c(1,5,0.1,2), mar=c(4,2,2,2))
for(s in 2:6){
  boxplot(bias ~ met, data = box_dat[box_dat$sce==s, ], las = 1, 
          ylab = "", xlab = "", xaxt = "n", ylim = c(0, 0.4),
          beside=TRUE,  cex.axis = 1.25,
          main = paste("S", s-1, sep=""), cex.lab = 1.5)
  axis(1, at = seq(1, 7, 1), labels = FALSE)
  labels <- c("REF", "Pool", "Test&Pool",  "BMA", "EXNEX", "BMA-2D", "EXNEX-2D")
  text(x = seq(1, 7, 1), y = rep(par("usr")[3] - 0.015, 7), srt = 45, adj = 1, 
       labels = labels, xpd = TRUE, cex = 1, font = 2)
  
}
dev.off()

bias_name_2 <- paste0("NI_bias_2_","r_",r,".pdf")
pdf(bias_name_2, width = 12, height = 8)
par(mfrow=c(2, 3), oma=c(1,5,0.1,2), mar=c(4,2,2,2))
for(s in 7:12){
  boxplot(bias ~ met, data = box_dat[box_dat$sce==s, ], las = 1, 
          ylab = "", xlab = "", xaxt = "n", ylim = c(0, 0.4),
          beside=TRUE,  cex.axis = 1.25, 
          main = paste("S", s-1, sep=""), cex.lab = 1.5)
  axis(1, at = seq(1, 7, 1), labels = FALSE)
  labels <- c("REF", "Pool", "Test&Pool",  "BMA", "EXNEX", "BMA-2D", "EXNEX-2D")
  text(x = seq(1, 7, 1), y = rep(par("usr")[3] - 0.015, 7), srt = 45, adj = 1, 
       labels = labels, xpd = TRUE, cex = 1, font = 2)
  
}
dev.off()
################################################################################















