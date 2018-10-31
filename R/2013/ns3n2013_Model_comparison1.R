###############################################################################################
#   This program plots legal abundance from different model configurations
###############################################################################################
library(PBSadmb)
library(gplots)
setwd('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2013/Model')
###############################################################################################
#   Crab Data construction 
###############################################################################################
rm(list=ls(all=TRUE))
admodeldir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2013/Model/')
setwd(admodeldir)

###############################################################################################
#   1.0 Read source file 
###############################################################################################
# Read read.R for program use
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/read_crab_data2013.R')
# Crab ouptutdata reading routine 
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/read_rep_data2.R')

##############################################################################################
#   1.1 Define input file name and output file name and locations. 
###############################################################################################
# infun: input file name 
datadir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/')
infn <- paste(datadir,'ns3n2012R.txt',sep='')
# outfn: output file name  
outfn <- paste(datadir,'ns3n2013R.dat',sep='')
pdffile <- paste(admodeldir,'ns3n2013Rcpue.pdf',sep='')


###############################################################################################
#   1.1.1 ADMB program file locations. 
###############################################################################################
# ADMB model file name
prefix1 <- paste(admodeldir,'Original/S3-1+qnhalf/','ns3n2013cpue',sep='')
prefix2 <- paste(admodeldir,'Original/S3-6+qnhalf+M0.24/','ns3n2013cpue',sep='')
prefix3 <- paste(admodeldir,'Original/S3-7+qnhalf+M0.30/','ns3n2013cpue',sep='')
prefix4 <- paste(admodeldir,'Original/S4-3-wocpue-nhalf-M0.30/','ns3n2013cpue',sep='')
prefix5 <- paste(admodeldir,'Original/S4-4-wocpue-q/','ns3n2013cpue',sep='')
prefix6 <- paste(admodeldir,'Original/S4-5-wocpue-qnhalf/','ns3n2013cpue',sep='')
prefix7 <- paste(admodeldir,'Original/S4-6-wocpue-qnhalf-M0.24/','ns3n2013cpue',sep='')
prefix8 <- paste(admodeldir,'Original/S4-7-wocpue-qnhalf-M0.30/','ns3n2013cpue',sep='')

legendtext <- c('S3-1','S3-6','S3-7')
# fn: ADMB report ouput file name 
# fnpar: ADMB parameter ouput file name 
fn <- paste(prefix,'.rep',sep='')
fnpar <- paste(prefix,'.std',sep='')

###############################################################################################
#   1.2 Red ADMB data into R  
###############################################################################################
# Read crab data into R List file 
dat<-datatoR(infn)
dat$maxss <- 20
dat$lamc <- 1
#dat$cv3[1:2] <- 0.01
#dat$ms <- 1
dat$M <- 0.30
dat$slt6 <- 1

########  Add weights to Summer pot survey likelihood ############################################
dat$log_initp  <- 8
dat$initp_phase <- 1
dat$wsp <- c(0,0)
dat$qtno_phase <- 1
dat$qtad_phase <- -1 


###############################################################################################
#   1.3 Read output Data from multiple sources 
###############################################################################################
rept1 <- reptoRlist(paste(prefix1,'.rep',sep=''))
rept2 <- reptoRlist(paste(prefix2,'.rep',sep=''))
rept3 <- reptoRlist(paste(prefix3,'.rep',sep=''))
rept4 <- reptoRlist(paste(prefix4,'.rep',sep=''))
rept5 <- reptoRlist(paste(prefix5,'.rep',sep=''))
rept6 <- reptoRlist(paste(prefix6,'.rep',sep=''))
rept7 <- reptoRlist(paste(prefix7,'.rep',sep=''))
rept8 <- reptoRlist(paste(prefix8,'.rep',sep=''))
par1 <- read.table(paste(prefix1,'.std',sep=''), sep="",header=T)
pars1 <- data.frame(t(par1[,3]))
names(pars1) <- t(par1[,2])
par2 <- read.table(paste(prefix2,'.std',sep=''), sep="",header=T)
pars2 <- data.frame(t(par2[,3]))
names(pars2) <- t(par2[,2])
par3 <- read.table(paste(prefix3,'.std',sep=''), sep="",header=T)
pars3 <- data.frame(t(par3[,3]))
names(pars3) <- t(par3[,2])
par4 <- read.table(paste(prefix4,'.std',sep=''), sep="",header=T)
pars4 <- data.frame(t(par4[,3]))
names(pars4) <- t(par4[,2])
par5 <- read.table(paste(prefix5,'.std',sep=''), sep="",header=T)
pars5 <- data.frame(t(par5[,3]))
names(pars5) <- t(par5[,2])
par6 <- read.table(paste(prefix6,'.std',sep=''), sep="",header=T)
pars6 <- data.frame(t(par6[,3]))
names(pars6) <- t(par6[,2])
par7 <- read.table(paste(prefix7,'.std',sep=''), sep="",header=T)
pars7 <- data.frame(t(par7[,3]))
names(pars7) <- t(par7[,2])
par8 <- read.table(paste(prefix8,'.std',sep=''), sep="",header=T)
pars8 <- data.frame(t(par8[,3]))
names(pars8) <- t(par8[,2])


###############################################################################################
#   1.4  plot 
###############################################################################################
windows(record=TRUE)
year <- seq(1976,1975+dat$ny)
# Figure 1. Total estimated abundance 
total.abundance1 <- colSums(rept1$nps+rept1$ops)/1000
total.abundance2 <- colSums(rept2$nps+rept2$ops)/1000
total.abundance3 <- colSums(rept3$nps+rept3$ops)/1000
total.abundance4 <- colSums(rept4$nps+rept1$ops)/1000
total.abundance5 <- colSums(rept5$nps+rept2$ops)/1000
total.abundance6 <- colSums(rept6$nps+rept3$ops)/1000
total.abundance7 <- colSums(rept7$nps+rept1$ops)/1000
total.abundance8 <- colSums(rept8$nps+rept2$ops)/1000

# Plot trawl survey data with CI
ciw <- 2*dat$tt*dat$cv3/1000
plotCI(year[dat$it],dat$tt/1000,uiw=ciw,ylab = 'Total Crab Abundance (million)', xlab = 'Year', ylim= c(0,8),bty='l')
lines(year,total.abundance1[1:dat$ny],lwd=2,col=1)
lines(year,total.abundance2[1:dat$ny],lwd=1,col=2)
lines(year,total.abundance3[1:dat$ny],lwd=1,col=3)
lines(year,total.abundance4[1:dat$ny],lwd=1,col=4)
lines(year,total.abundance5[1:dat$ny],lwd=1,col=5)
lines(year,total.abundance6[1:dat$ny],lwd=1,col=6)
lines(year,total.abundance7[1:dat$ny],lwd=1,col=7)
lines(year,total.abundance8[1:dat$ny],lwd=1,col=8)
#points(year[dat$it],etrawl/1000)
legend("topright", lwd = c(2,1,1,1,1,1,1,1), col=c(1,2,3,4,5,6,7,8), legend = legendtext)
title(main='Trawl survey crab abundance')

# Figure 2.1 Total, Recruites, Legals estimated abundance 
plot(year, rept1$Legal/1000,ylab = 'Abundance (million crabs)', typ = 'l', lwd=2, xlab = 'Year',ylim= c(0,6),bty='l')
lines(year, rept2$Legal/1000,lwd=2,col=2, lty=2)
lines(year, rept3$Legal/1000,lwd=2,col=4, lty=4)
legend("topright", lty = c(1,2,4), lwd = c(2,2,2), col=c(1,2,4), legend = legendtext)
title(main='Trawl survey Legal crab abundance')

# Figure 2.2 Total, Recruites, Legals estimated abundance 
plot(year[1:(dat$ny-1)], rept1$Rec[1:(dat$ny-1)]/1000,ylab = 'Abundance (million crabs)', typ = 'l', lwd=2, xlab = 'Year',ylim= c(0,6),bty='l')
lines(year[1:(dat$ny-1)], rept2$Rec[1:(dat$ny-1)]/1000,lwd=2,col=2,lty=2)
lines(year[1:(dat$ny-1)], rept3$Rec[1:(dat$ny-1)]/1000,lwd=2,col=4,lty=4)
legend("topright", lty = c(1,2,4), lwd = c(2,2,2), col=c(1,2,4), legend = legendtext)
title(main='Trawl survey Recruit crab abundance')

# Figure 3. cpue data 
ciw <- 2*dat$stcpuese[-dat$ny]
plotCI(year[-dat$ny], dat$stcpue[-dat$ny], uiw=ciw, ylab = 'Standardized CPUE',  xlab = 'Year', ylim=c(0,7), gap=0.2, bty='l')
lines(year[-dat$ny],rept1$ecpue[-dat$ny],lwd=2)
lines(year[-dat$ny],rept2$ecpue[-dat$ny],lwd=2, col=2, lty=2)
lines(year[-dat$ny],rept3$ecpue[-dat$ny],lwd=2, col=4, lty=4)
legend("topright",lty =c(0,1,2,4), lwd = c(0,2,2,2), pch=c(1,NA,NA,NA), col=c(1,1,2,4), legend = c('Observed','S3-1','S3-6','S3-7'))
title(main='Summer commercial standardized cpue')


# Figure 4. Catch data 
par(mar=c(5,4,4,5))
tcatch <- dat$tc[-dat$ny]+dat$twc+dat$tsc
hrate1 <- tcatch/rept1$Legal[-dat$ny]
hrate2 <- tcatch/rept2$Legal[-dat$ny]
hrate3 <- tcatch/rept3$Legal[-dat$ny]
plot(year[-dat$ny], tcatch[-dat$ny]/1000, ylab = 'Total Catch (million)', lwd=2, typ = 'l', xlab = 'Year', bty='u')
par(new=TRUE)
plot(year[-dat$ny],hrate1,lty=5, lwd=2, typ = 'l',bty='l', xaxt='n', yaxt='n', xlab='',ylab='', ylim = c(0,max(hrate1)))
lines(year[-dat$ny],hrate2,lty=2, lwd=2, col=2)
lines(year[-dat$ny],hrate3,lty=4, lwd=2, col=4)
axis(4)
mtext("Estimated harvest rate",side=4,line=3)
legend("topright", lty = c(1,2), lwd = 2, legend=c('Total Catch', 'Estimated Harvest Rate'))
title(main='Total catch & Estimated harvest rate')

# Figure 5. Catch data 
mmb1 <- rept1$MMB
mmb2 <- rept2$MMB
mmb3 <- rept3$MMB
par(mfrow=c(3,3),mar = c(2, 2, 2, 2))
plot(mmb1[1:(dat$ny-1)]/1000, hrate1, main='S3-1', ylab = 'Harvest rate', typ = 'p', xlab = 'MMB',xlim=c(0,16),ylim=c(0,0.5))
plot(mmb2[1:(dat$ny-1)]/1000, hrate2, main='S3-6',  col =2,xlim=c(0,16),ylim=c(0,0.5))
plot(mmb3[1:(dat$ny-1)]/1000, hrate3, main='S3-7',  col =4,xlim=c(0,16),ylim=c(0,0.5))

plot(mmb2[1:(dat$ny-1)]/1000, hrate2,ylab = 'Harvest rate', typ = 'p',pch = 20, xlab = 'MMB')
plot(mmb3[1:(dat$ny-1)]/1000, hrate3,ylab = 'Harvest rate', typ = 'p',pch = 20, xlab = 'MMB')





##############################################################################################################
#  Molting probability and Net Selectivity
##############################################################################################################
par(mfrow=c(2,3),mar = c(2, 1.5, 1.5, 1.5),oma=c(4,1,2,1))
sizeclass <- c(78.5, 88.5, 98.5, 108.5, 118.5, 128.5, 128.5)
lengthclass  <- c(74.0, 84.0, 94.0, 104, 114, 124, 134)
moltp <- 1-1/(1+exp(-(exp(pars1$log_mo1)*(sizeclass-exp(pars1$log_mo2)))))
moltps1 <- moltp/moltp[1]
moltp <- 1-1/(1+exp(-(exp(pars2$log_mo1)*(sizeclass-exp(pars2$log_mo2)))))
moltps2 <- moltp/moltp[1]
moltp <- 1-1/(1+exp(-(exp(pars3$log_mo1)*(sizeclass-exp(pars3$log_mo2)))))
moltps3 <- moltp/moltp[1]
plot(lengthclass, moltps1,type='s', main='Molting Probability', ylim=c(0,1))
lines(lengthclass, moltps2,type='s', col=2, lty=2, ylim=c(0,1))
lines(lengthclass, moltps3,type='s', col=4, lty=4, ylim=c(0,1))

trawlp <- 1/(1+exp(-(exp(pars1$log_st1)*(sizeclass-exp(pars1$log_st2)))))
trawlps1 <- ifelse(trawlp/trawlp[5]>1,1,trawlp/trawlp[5])
trawlp <- 1/(1+exp(-(exp(pars2$log_st1)*(sizeclass-exp(pars2$log_st2)))))
trawlps2 <- ifelse(trawlp/trawlp[5]>1,1,trawlp/trawlp[5])
trawlp <- 1/(1+exp(-(exp(pars3$log_st1)*(sizeclass-exp(pars3$log_st2)))))
trawlps3 <- ifelse(trawlp/trawlp[5]>1,1,trawlp/trawlp[5])
plot(lengthclass, trawlps1,type='s', main='Trawl Selectivity', ylim=c(0,1))
lines(lengthclass, trawlps2,type='s', col=2, lty=2, ylim=c(0,1))
lines(lengthclass, trawlps3,type='s', col=4, lty=4, ylim=c(0,1))

wpotp <- 1/(1+exp(-(exp(pars1$log_sw1)*(sizeclass-exp(pars1$log_sw2)))))
wpotps1 <- ifelse(wpotp/wpotp[5]>1,1,wpotp/wpotp[5])
wpotps1[6:7] <- (pars1$sw3)
wpotp <- 1/(1+exp(-(exp(pars2$log_sw1)*(sizeclass-exp(pars2$log_sw2)))))
wpotps2 <- ifelse(wpotp/wpotp[5]>1,1,wpotp/wpotp[5])
wpotps2[6:7] <- (pars2$sw3)
wpotp <- 1/(1+exp(-(exp(pars3$log_sw1)*(sizeclass-exp(pars3$log_sw2)))))
wpotps3 <- ifelse(wpotp/wpotp[5]>1,1,wpotp/wpotp[5])
wpotps3[6:7] <- (pars3$sw3)
plot(lengthclass, wpotps1,type='s', main='Winter pot Selectivity', ylim=c(0,1))
lines(lengthclass, wpotps2,type='s', col=2, lty=2, ylim=c(0,1))
lines(lengthclass, wpotps3,type='s', col=4, lty=4, ylim=c(0,1))

sc1p <- 1/(1+exp(-(exp(pars1$log_sc1)*(sizeclass-exp(pars1$log_sc2)))))
sc1ps1 <- ifelse(sc1p/sc1p[5]>1,1,sc1p/sc1p[5])
sc1p <- 1/(1+exp(-(exp(pars2$log_sc1)*(sizeclass-exp(pars2$log_sc2)))))
sc1ps2 <- ifelse(sc1p/sc1p[5]>1,1,sc1p/sc1p[5])
sc1p <- 1/(1+exp(-(exp(pars3$log_sc1)*(sizeclass-exp(pars3$log_sc2)))))
sc1ps3 <- ifelse(sc1p/sc1p[5]>1,1,sc1p/sc1p[5])
plot(lengthclass, sc1ps1,type='s', main='Commercial 77-92 Selectivity', ylim=c(0,1))
lines(lengthclass, sc1ps2,type='s', col=2, lty=2, ylim=c(0,1))
lines(lengthclass, sc1ps3,type='s', col=4, lty=4, ylim=c(0,1))

sc2p <- 1/(1+exp(-(exp(pars1$log_sc3)*(sizeclass-exp(pars1$log_sc4)))))
sc2ps1 <- ifelse(sc2p/sc2p[5]>1,1,sc2p/sc2p[5])
sc2p <- 1/(1+exp(-(exp(pars2$log_sc3)*(sizeclass-exp(pars2$log_sc4)))))
sc2ps2 <- ifelse(sc2p/sc2p[5]>1,1,sc2p/sc2p[5])
sc2p <- 1/(1+exp(-(exp(pars3$log_sc3)*(sizeclass-exp(pars3$log_sc4)))))
sc2ps3 <- ifelse(sc2p/sc2p[5]>1,1,sc2p/sc2p[5])
plot(lengthclass, sc2ps3,type='s', main='Commercial 93-12 Selectivity', ylim=c(0,1))
lines(lengthclass, sc2ps2,type='s', col=2, lty=2, ylim=c(0,1))
lines(lengthclass, sc2ps3,type='s', col=4, lty=4, ylim=c(0,1))


##############################################################################################################
#  Figure for Length Frequency Plot
##############################################################################################################
par(mfrow=c(2,3),mar = c(2, 1.5, 1.5, 1.5),oma=c(4,1,2,1))
lengthclass  <- c(74.0, 84.0, 94.0, 104, 114, 124,134)
# Create commercial length vs predicted length frequencies
cdata <- rowSums((dat$onc+dat$ooc)*dat$tc)
ccdata <- c(cumsum(cdata)/sum(cdata),1)
pdata <- rowSums((rept1$enc + rept1$eoc)*dat$tc)
cpdata1 <- c(cumsum(pdata)/sum(pdata),1)
pdata <- rowSums((rept2$enc + rept2$eoc)*dat$tc)
cpdata2 <- c(cumsum(pdata)/sum(pdata),1)
pdata <- rowSums((rept3$enc + rept3$eoc)*dat$tc)
cpdata3 <- c(cumsum(pdata)/sum(pdata),1)
plot(lengthclass, ccdata, type='s', main='Commercial Harvest', xlab = 'Length', ylab='')
lines(lengthclass, cpdata1, type='s',lty=5, col=1)
lines(lengthclass, cpdata2, type='s',lty=2, col=2)
lines(lengthclass, cpdata3, type='s',lty=4, col=4)

# Create Winter length vs predicted length freq figure 
cdata <- rowSums((dat$onw+dat$oow)*dat$sw)
ccdata <- c(cumsum(cdata)/sum(cdata),1)
pdata <- rowSums((rept1$enw + rept1$eow)*dat$sw)
cpdata1 <- c(cumsum(pdata)/sum(pdata),1)
pdata <- rowSums((rept2$enw + rept2$eow)*dat$sw)
cpdata2 <- c(cumsum(pdata)/sum(pdata),1)
pdata <- rowSums((rept3$enw + rept3$eow)*dat$sw)
cpdata3 <- c(cumsum(pdata)/sum(pdata),1)
plot(lengthclass, ccdata, type='s', main='Winter Pot Survey', xlab = 'Length', ylab='')
lines(lengthclass, cpdata1, type='s',lty=5, col=1)
lines(lengthclass, cpdata2, type='s',lty=2, col=2)
lines(lengthclass, cpdata3, type='s',lty=4, col=4)

# Create Trawl length vs predicted length freq figure 
cdata <- rowSums((dat$ont+dat$oot)*dat$tt)
ccdata <- c(cumsum(cdata)/sum(cdata),1)
pdata <- rowSums((rept1$ent)*rept1$ett)
cpdata1 <- c(cumsum(pdata)/sum(pdata),1)
pdata <- rowSums((rept2$ent)*rept2$ett)
cpdata2 <- c(cumsum(pdata)/sum(pdata),1)
pdata <- rowSums((rept3$ent)*rept3$ett)
cpdata3 <- c(cumsum(pdata)/sum(pdata),1)
plot(lengthclass, ccdata, type='s', main='Trawl Survey', xlab = 'Length', ylab='')
lines(lengthclass, cpdata1, type='s',lty=5, col=1)
lines(lengthclass, cpdata2, type='s',lty=2, col=2)
lines(lengthclass, cpdata3, type='s',lty=4, col=4)

# Create Observer survey length vs predicted length freq figure 
cdata <- rowSums((dat$ono+dat$ooo)*dat$so)
ccdata <- c(cumsum(cdata)/sum(cdata),1)
pdata <- rowSums((rept1$eno + rept1$eoo)*dat$so)
cpdata1 <- c(cumsum(pdata)/sum(pdata),1)
pdata <- rowSums((rept2$eno + rept2$eoo)*dat$so)
cpdata2 <- c(cumsum(pdata)/sum(pdata),1)
pdata <- rowSums((rept3$eno + rept3$eoo)*dat$so)
cpdata3 <- c(cumsum(pdata)/sum(pdata),1)
plot(lengthclass, ccdata, type='s', main='Summer Observer', xlab = 'Length', ylab='')
lines(lengthclass, cpdata1, type='s',lty=5, col=1)
lines(lengthclass, cpdata2, type='s',lty=2, col=2)
lines(lengthclass, cpdata3, type='s',lty=4, col=4)




# Create Pot survey length vs predicted length freq figure 
cdata <- rowSums((dat$onp+dat$oop)*dat$tp)
ccdata <- c(cumsum(cdata)/sum(cdata),1)
pdata <- rowSums((rept$enp + rept$eop)*rept$etp)
cpdata <- c(cumsum(pdata)/sum(pdata),1)
plot(lengthclass, ccdata, type='s', main='Summer Pot Survey', xlab = 'Length', ylab='')
lines(lengthclass, cpdata, type='s', lty=5, col='red')
mtext("Length frequency: Observed (black) vs. Predicted (red)", side=3,outer=T)


