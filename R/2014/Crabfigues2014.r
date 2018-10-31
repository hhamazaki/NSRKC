##############################################################################################################
#  Figure producing R function rutine
##############################################################################################################
library(gplots)
library(PBSadmb)
rm(list=ls(all=TRUE))
admodeldir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2014/Model/')
setwd(admodeldir)

###############################################################################################
#   1.0 Read source file 
###############################################################################################
# Read read.R for program use
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/read_crab_data2013.R')
# Crab ouptutdata reading routine 
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/read_rep_data2.R')
# Crab Julyschedule change routine
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/make_July_data.R')

##############################################################################################
#   1.1 Define input file name and output file name and locations. 
###############################################################################################
# infun: input file name 
datadir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/')
infn <- paste(datadir,'ns3n2014R.txt',sep='')
# outfn: output file name  

##############################################################################################
#   Activate either, based on model you use 
###############################################################################################
outfn <- paste(datadir,'ns3n2014R.dat',sep='')
# July model
#outfn <- paste(datadir,'ns3n2014RJ.dat',sep='')

pdffile <- paste(admodeldir,'ns3n2014.pdf',sep='')

###############################################################################################
#   1.1.1 ADMB program file locations: Choose prefix based on your modeling 
###############################################################################################
# ADMB model file name
prefix <- paste(admodeldir,'ns3n2014_1_4c',sep='')
# July model 
#prefix <- paste(admodeldir,'ns3n2014_Jul1',sep='')
# fn: ADMB report ouput file name 

# fnpar: ADMB parameter ouput file name 
fn <- paste(prefix,'.rep',sep='')
fnpar <- paste(prefix,'.std',sep='')


###############################################################################################
#   1.2 Red ADMB data into R  
###############################################################################################
# Read crab data into R List file 
dat<-datatoR(infn)
# Convert to July model data 
#dat<-Julydata(dat)

#dat$cv3[1:2] <- 0.01
#dat$ms <- 1.0
#dat$M <- 0.34
#dat$slt6 <- 1
#dat$ny <- 37
#dat$maxss <- 2
#dat$maxsc <- 1
########  Add weights to Summer pot survey likelihood ############################################
dat$SDRec <- 0.5
dat$ms_phase <- -1
dat$lamc <- 1
dat$lamw <- 1


########  Remove last observer data ############################################
#dat$ono <- dat$ono[,-dat$nyo]
#dat$ooo <- dat$ooo[,-dat$nyo]
#dat$io <- dat$io[-dat$nyo]
#dat$so <- dat$so[-dat$nyo]
#dat$nyo <- dat$nyo-1

########  Remove cpue data ############################################
#dat$stcpue[dat$ny] <- 0 
#dat$stcpuese[dat$ny] <- 0 

  
######### Add 0 to projected year ##################################################
#dat$ny <- dat$ny+1
#dat$te[dat$ny] <- 0
#dat$stcpue[dat$ny] <- 0
#dat$stcpuese[dat$ny] <- 0
#dat$sc[dat$ny] <- 0
#dat$tc[dat$ny] <- 0
#dat$ys[dat$ny] <- 0
#dat$twc[dat$ny] <- 0
#dat$tsc[dat$ny] <- 0
#dat$tsct[dat$ny] <- 0
#dat$onc <- cbind(dat$onc,c(0,0,0,0,0,0))
#dat$ooc <- cbind(dat$ooc,c(0,0,0,0,0,0))




###################################################################################################	
	
	
makedata <- function(dat,outfn){
# Save original file as filename_0.dat
file.remove(paste(outfn))
out_file <- file(paste(outfn), open='a')
for (j in seq_along(dat)){
    write.table(paste("#",names(dat)[j]), file=out_file, sep="", dec=".", quote=FALSE, col.names=FALSE, row.names=FALSE)  
	write(t(dat[[j]]), file=out_file, sep=" ", 
	ncolumns= ifelse(is.numeric(dim(dat[[j]])[2]),dim(dat[[j]])[2],length(dat[[j]])))
	}
close(out_file) 
}
makedata(dat,outfn)

# Run the ADMB 
argvec <- paste(' -ind ',outfn,sep="")
runAD(prefix=prefix, argvec = argvec, logfile=TRUE, add = FALSE, verbose=TRUE)

# Read output file
rept <- reptoRlist(fn)
#rept
par1 <- read.table(fnpar, sep="",header=T)
pars <- data.frame(t(par1[,3]))
names(pars) <- t(par1[,2])
rept$f
rept$tf
par1
test <- paste(prefix,'.par',sep="") 
test1 <- as.double(scan(test,nlines=1,quiet=T,what="")) 
maxgrad <-test1[16]
maxgrad

##############################################################################################################
#  Calculate RMSE and Residual plots
##############################################################################################################
trawl1 <- ifelse(is.null(pars$qtno),1,pars$qtno)*(rept$ett1[1:6]+rept$ett2[1:6])
trawl2 <- ifelse(is.null(pars$qtad),1,pars$qtad)*(rept$ett1[7:dat$nyt]+rept$ett2[7:dat$nyt])
etrawl <- c(trawl1,trawl2)
liket <- (log(dat$tt+0.001)-log(etrawl+0.001))^2/log(dat$cv3*dat$cv3+1)
RMSE.trawl <- sqrt(sum((log(dat$tt)-log(etrawl))^2)/length(dat$tt))
RMSE.pot <- sqrt(sum((log(dat$tp)-log(rept$etp1+rept$etp2))^2)/length(dat$tp))
RMSE.cpue <- sqrt(sum((log(dat$stcpue[-c(1,13,16,dat$ny)])-log(rept$ecpue[-c(1,13,16,dat$ny)]))^2)/length(dat$stcpue[-c(1,13, 16,dat$ny)]))
RMSE.trawl
RMSE.pot
RMSE.cpue
resd.trawl <- dat$tt - etrawl
resd.pot <- dat$tp - (rept$etp1+rept$etp2)
resd.cpue <- dat$stcpue[-c(1,13,16,dat$ny)] - rept$ecpue[-c(1,13,16,dat$ny)]

##############################################################################################################
#  Figures 
##############################################################################################################

windows(record=TRUE)# Create year data 
year <- seq(1976,1975+dat$ny+1)
lenclass <- c('74-83','84-93','94-103','104-113','114-123','>123')

#dev.copy(pdf, pdffile)
#dev.copy(win.metafile, metafile)
##############################################################################################################
#  Molting probability and Net Selectivity
##############################################################################################################
par(mfrow=c(2,3),mar = c(2, 1.5, 1.5, 1.5),oma=c(4,1,2,1))
sizeclass <- c(78.5, 88.5, 98.5, 108.5, 118.5, 128.5, 128.5)
lengthclass  <- c(74.0, 84.0, 94.0, 104, 114, 124, 134)
moltp <- 1-1/(1+exp(exp(pars$log_mo1)*(78.5-sizeclass)+log(1.0/0.01-1.0)))
plot(lengthclass, moltp,type='s', main='Molting Probability', ylim=c(0,1))
trawlpn <- 1/(1+exp(exp(pars$log_st1)*(128.5-sizeclass)+log(1.0/0.99-1.0)))
plot(lengthclass, trawlpn,type='s', main='NOAA Trawl Selectivity', ylim=c(0,1))
trawlpa <- 1/(1+exp(exp(pars$log_st3)*(128.5-sizeclass)+log(1.0/0.99-1.0)))
plot(lengthclass, trawlpa,type='s', main='ADFG Trawl Selectivity', ylim=c(0,1))
#spotp <- 1/(1+exp(-(exp(pars$log_sp1)*(sizeclass-exp(pars$log_sp2)))))
#spotps <- ifelse(spotp/spotp[5]>1,1,spotp/spotp[5])
#plot(lengthclass, spotps,type='s', main='Summer pot Selectivity', ylim=c(0,1))
wpotp <- 1/(1+exp(exp(pars$log_sw1)*(128.5-sizeclass)+log(1.0/0.99-1.0)))
wpotp[6:7] <- (pars$sw3)
plot(lengthclass, wpotp,type='s', main='Winter pot Selectivity', ylim=c(0,1))
sc1p <- 1/(1+exp(exp(pars$log_sc1)*(128.5-sizeclass)+log(1.0/0.99-1.0)))
plot(lengthclass, sc1p,type='s', main='Commercial 77-92 Selectivity', ylim=c(0,1))
#sc2p <- 1/(1+exp(exp(pars$log_sc3)*(128.5-sizeclass)+log(1.0/0.99-1.0)))
#plot(lengthclass, sc2p,type='s', main='Ctommercial 93-05 Selectivity', ylim=c(0,1))
sc3p <- 1/(1+exp(exp(pars$log_sc5)*(128.5-sizeclass)+log(1.0/0.99-1.0)))
plot(lengthclass, sc3p,type='s', main='Commercial 93-12 Selectivity', ylim=c(0,1))


##############################################################################################################
#  Plot residual QQ plot for survey
##############################################################################################################
par(mfrow=c(3,3),mar = c(2, 2, 2, 2),oma=c(3,1,2,1))
hist(resd.trawl, main='Trawl survey')
qqnorm(resd.trawl)
qqline(resd.trawl)
plot((rept$ett1+rept$ett2),resd.trawl,main='predicted vs. residual')
abline(h=0)
hist(resd.cpue,, main='Commercial CPUE')
qqnorm(resd.cpue)
qqline(resd.cpue)
plot(rept$ecpue[-c(1,13,16,dat$ny)],resd.cpue,main='predicted vs. residual')
abline(h=0)
#hist(resd.pot,main='Summer pot survey Residuals')
#qqnorm(resd.pot)
#qqline(resd.pot)
#plot((rept$etp1+rept$etp2),resd.pot)
#abline(h=0)
mtext("Residuals Histogram, Q-Q Plot, Predicted vs. Residual", side=3, outer=T)

##############################################################################################################
#  Calculate and Plot Effective sample size 
##############################################################################################################
#trawl survey effective sample size
ot <- dat$ont+dat$oot
tefn <-colSums(rept$ent*(1-rept$ent))/colSums((ot-rept$ent)^2)
tefn
tefnn <- ifelse(dat$st*0.5>dat$maxss,dat$maxss,dat$st*0.5)
 
#Summer pot survey effective sample size
#spp <- rept$enp+rept$eop
#sp <- dat$onp+dat$oop
#spefn <-colSums(spp*(1-spp))/colSums((spp-sp)^2)
#spefnn <- ifelse(dat$sp*0.5>dat$maxss,dat$maxss,dat$sp*0.5)

#commercial fishery effective sample size
scp <- rept$enc+rept$eoc
sc <- dat$onc+dat$ooc
scefn <-colSums(scp*(1-scp))/colSums((scp-sc)^2)
scefn <-scefn[-c(1,16)]
scefnn <- ifelse(dat$sc*0.1>dat$maxss/2,dat$maxss/2,dat$sc*0.1)[-c(1,16)]

#Winter Pot survey effective sample size
wpp <- rept$enw+rept$eow
wp <- dat$onw+dat$oow
wpefn <-colSums(wpp*(1-wpp))/colSums((wpp-wp)^2)
wpefnn <- ifelse(dat$sw*0.1>dat$maxss/2,dat$maxss/2,dat$sw*0.1)

# Observer survey effective sample size
opp <- rept$eno+rept$eoo
op <- dat$ono+dat$ooo
opefn <-colSums(opp*(1-opp))/colSums((opp-op)^2)
opefnn <- ifelse(dat$so*0.1>dat$maxss/2,dat$maxss/2,dat$so*0.1)
par(mfrow=c(5,3),mar = c(2, 2, 2, 2),oma=c(3,1,2,1))

# Trawl Survey Figure
hist(tefn,breaks=10, main='Trawl survey')
abline(v=c(dat$maxss),lwd=2)
plot(tefnn,tefn, main='Trawl survey', xlim = c(0,dat$maxss),ylim=c(0,max(tefn)))
abline(lm(tefn~ tefnn -1),lty=2)
lines(c(0,dat$maxss),c(0,dat$maxss))
plot(dat$it+1975, tefn)

# Commercial Catch Figure
hist(scefn,breaks=10, main='Commercial Catch')
abline(v=c(dat$maxss/2),lwd=2)
plot(scefnn, scefn, main='Commercial Catch', xlim = c(0,dat$maxss/2),ylim=c(0,max(scefn)))
lines(scefn~scefnn-1,lty=2)
lines(c(0,dat$maxss/2),c(0,dat$maxss/2))
plot(year[-c(1,16,dat$ny+1)], scefn)

# Winter Pot Survey Figure
hist(wpefn,breaks=5, main='Winter pot survey')
abline(v=c(dat$maxss/2),lwd=2)
plot(wpefnn, wpefn, main='Winter pot survey', xlim = c(0,dat$maxss/2),ylim=c(0,max(wpefn)))
abline(lm(wpefn~wpefnn-1),lty=2)
lines(c(0,dat$maxss/2),c(0,dat$maxss/2))
plot(dat$iw+1975, wpefn)

# Observer Survey Figure
hist(opefn,breaks=10, main='Observer survey')
abline(v=c(dat$maxss/2),lwd=2)
plot(opefnn, opefn, main='Observer survey', xlim = c(0,dat$maxss/2),ylim=c(0,max(opefn)))
abline(lm(opefn~opefnn-1),lty=2)
lines(c(0,dat$maxss/2),c(0,dat$maxss/2))
plot(dat$io+1975, opefn)

# Summer Pot Survey Figure
#hist(spefn,breaks=10, main='Summer pot survey')
#abline(v=c(dat$maxss),lwd=2)
#plot(spefnn,spefn, main='Summer pot survey', xlim = c(0,dat$maxss),ylim=c(0,max(spefn)))
#abline(lm(spefn ~ spefnn -1),lty=2)
#lines(c(0,dat$maxss),c(0,dat$maxss))
#plot(dat$ip+1975, spefn)
mtext("Effective sample size", side=3, outer=T)

par(mfrow=c(1,1), mar=c(5,5,4,4))
##############################################################################################################
#  Figure for trends.
##############################################################################################################

# Figure 1. Total estimated abundance 
# Plot trawl survey data with CI
ciw <- 2*dat$tt*dat$cv3/1000
plotCI(year[dat$it],dat$tt/1000,uiw=ciw,ylab = 'Crab Abundance (million)', xlab = 'Year', ylim= c(0,8),bty='l')
#points(year[dat$ip], dat$tp/1000, pch= 19)
# Model estiamted trawl survey abundance
points(year[dat$it],etrawl/1000,pch=19, col=4)
#lines(year,total.abundance[1:dat$ny],lwd=2)
legend("topright", pch = c(1,20,20), col= c(1,4,1), legend = c('Observed','Predicted'))
title(main='Trawl survey crab abundance')


# Figure 2. Total, Recruites, Legals estimated abundance 
#: Add projected number 
#rept$nps[1:5,dat$ny+1] <- rept$npw[1:5,dat$ny+1]*exp(-0.417*dat$M)
#rept$nps[6,dat$ny+1] <- rept$npw[6,dat$ny+1]*exp(-0.417*dat$M*dat$ms)
#rept$ops[1:5,dat$ny+1] <- rept$opw[1:5,dat$ny+1]*exp(-0.417*dat$M)
#rept$ops[6,dat$ny+1] <- rept$opw[6,dat$ny+1]*exp(-0.417*dat$M*dat$ms)

total.abundance <- colSums(rept$nps+rept$ops)/1000
plot(year[1:(dat$ny+1)], total.abundance[1:(dat$ny+1)], ylab = 'Abundance (million crabs)', typ = 'l', lwd=2, xlab = 'Year',ylim= c(0,10),bty='l')
legal.abundance <- colSums(rept$nps[4:6,]+rept$ops[4:6,])
lines(year[1:(dat$ny+1)], legal.abundance[1:(dat$ny+1)]/1000, lty=2)
recruit.abundance <- colSums(rept$nps[1:2,]+rept$ops[1:2,])
lines(year[1:(dat$ny+1)], recruit.abundance[1:(dat$ny+1)]/1000, lty=3)
legend("topright",lty = c(1,2,3), legend = c('total','legal','recruits'))
title(main='Modeled crab abundance July 01')


# Figure 5. MMB data 
mmb <- rept$MMB[1:(dat$ny+1)]
plot(year, mmb/1000, ylab = 'MMB (million lb)', lwd=2, typ = 'l', xlab = 'Year', bty='l')
bmsy <- mean(mmb[which(year==1980):(dat$ny+1)])
abline(h=bmsy/1000, lwd=2, lty=2)
cofl <- pars$last_mmb/bmsy
FOFL <- ifelse(cofl>1,dat$M,ifelse(0.25<cofl&&cofl<=1,dat$M*(cofl-0.1)/0.9,0))
OFLrb <- (1-exp(-FOFL))*pars$last_legal/1000
ABC <- OFLrb*0.9
title(main='MMB July 01')
tex <- c(paste('BMSY ',round(bmsy/1000,3),' mil.lb'),paste('MMB ',round(pars$last_mmb/1000,3),' mil.lb'),paste('Legal B ',round(pars$last_legal/1000,3),' mil.lb'),paste('FOFL ',round(FOFL,3)), paste('OFL ',round(OFLrb,3),' mil.lb'),paste('ABC ',round(ABC,3),' mil.lb'))
legend("topright",tex,bty='n',xjust=0) 


# Figure 3. Effort data 
#plot(year[-dat$ny], rept$ete[-dat$ny]/1000, ylab = 'Effort (1000 pot lifts)', lwd=2, typ = 'l', xlab = 'Year', ylim=c(0,40), bty='l')
#points(year[-dat$ny], dat$te[-dat$ny]/1000, pch=19)
#legend("topright",lty =c(1,0), lwd = c(2,0), pch=c(NA,19), legend = c('Predicted','Observed'))
#title(main='Summer commercial catch effort')
#sigmae <- sqrt(sum((log(dat$te[-c(1,16,dat$ny)])-log(rept$ete[-c(1,16,dat$ny)]))^2)/length(dat$te[-c(1,16,dat$ny)]))

# Figure 3. cpue data 
dat$stcpue[dat$stcpue==0]<-NA
ciw <- 2*dat$stcpuese
ciw2 <- 2*dat$stcpue*sqrt(exp(log((dat$stcpuese/dat$stcpue)^2+1)+pars$advar)-1)
ciw[ciw==0]<-NA
plot(year[-(dat$ny+1)], dat$stcpue, ylab = 'CPUE',  xlab = 'Year', ylim=c(0,7), bty='l')
arrows(year[-(dat$ny+1)],y0=dat$stcpue+ciw2,y1=dat$stcpue-ciw2,code=0,col='red')
arrows(year[-(dat$ny+1)],y0=dat$stcpue+ciw,y1=dat$stcpue-ciw,code=3,length=0.025, angle =90)
points(year[-(dat$ny+1)], dat$stcpue, pch=21, bg='white')
rept$ecpue[rept$ecpue==0]<-NA
lines(year[-(dat$ny+1)],rept$ecpue,lwd=2)
legend("topright",lty =c(1,0), lwd = c(2,0), pch=c(NA,19), legend = c('Predicted','Observed'))
title(main='Summer commercial standardized cpue')

# Figure 3.5. Length class data
par(mfrow=c(1,1), mar=c(5,5,4,4))
bdata <- rept$npw+rept$opw 
tdata <- t(dat$tt*t(ot))
plot(year, bdata[1,1:(dat$ny+1)], type='l', main='Estimated abundance by length class', ylab='Abundance', ylim=c(0,max(bdata)))
lines(year, bdata[2,1:(dat$ny+1)],lty=2)
lines(year, bdata[3,1:(dat$ny+1)],lty=1,col=2)
lines(year, bdata[4,1:(dat$ny+1)],lty=2,col=2)
lines(year, bdata[5,1:(dat$ny+1)],lty=1, col=4)
lines(year, bdata[6,1:(dat$ny+1)],lty=2, col=4)
legend("topright", bty ='n', lwd=1, lty=c(1,2,1,2,1,2),col=c(1,1,2,2,4,4),
legend = c('74-83','84-93','94-103','104-113','114-123','> 124'))
 

# Figure 4. Catch data 
par(mar=c(5,4,4,5))
tcatch <- dat$tc+dat$twc+dat$tws
hrate <- tcatch/rept$Legal[1:dat$ny]
grate <- c(0,(total.abundance[2:(dat$ny-1)]-total.abundance[1:(dat$ny-2)])/total.abundance[1:(dat$ny-2)])
plot(year[1:dat$ny], tcatch/1000, ylab = 'Total Catch (million)', lwd=2, typ = 'l', xlab = 'Year', bty='u')
par(new=TRUE)
plot(year[1:dat$ny],hrate[1:dat$ny],lty=2, lwd=2, typ = 'l',bty='l', xaxt='n', yaxt='n', xlab='',ylab='', ylim = c(0,max(hrate)))
axis(4)
mtext("Estimated harvest rate",side=4,line=3)
legend("topright", lty = c(1,2), lwd = 2, legend=c('Total Catch', 'Estimated Harvest Rate'))
title(main='Total catch & Harvest rate')


##############################################################################################################
#  Figure for size proportion 
##############################################################################################################
# Create commercial length vs predicted length freq figure 
par(mfrow=c(6,6),mar = c(1.2, 1.2, 1.2, 1.2),oma=c(3,1,2,1))
for (i in 2:15){
x <- plot(dat$onc[,i]+dat$ooc[,i],ylim=c(0,0.7),pch=19, cex=0.8, tck = -0.05, bty='l')
lines(rept$enc[,i]+rept$eoc[,i],lty=2)
text(2,0.6,paste(year[i]),cex=1.0)
}
for (i in 17:(dat$ny)){
x <- plot(dat$onc[,i]+dat$ooc[,i],ylim=c(0,0.7),pch=19, cex=0.8, tck = -0.05, bty='l')
lines(rept$enc[,i]+rept$eoc[,i],lty=2)
text(2,0.6,paste(year[i]),cex=1.0)
}
mtext("commercial harvest length: observed vs predicted", side=3, outer=T)
mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)

# Create winter length vs predicted length freq figure 
par(mfrow=c(6,6),mar = c(1.2, 1.2, 1.2, 1.2),oma=c(3,1,2,1))
for (i in 1:dat$nyw){
x <- plot(dat$onw[,i]+dat$oow[,i],ylim=c(0,0.5), pch=19, cex=0.8, tck = -0.05,bty='l', xlab = " ", ylab=" ")
lines(rept$enw[,i]+rept$eow[,i],lty=2)
text(2,0.45,paste(year[dat$iw[i]]),cex=1.0)
}
mtext("Winter pot length: observed vs predicted", side=3, outer=T)
mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)

# Create trawl length vs predicted length freq figure 
par(mfrow=c(6,6),mar = c(1.2, 1.2, 1.2, 1.2),oma=c(3,1,2,1))
for (i in 1:dat$nyt){
x <- plot(dat$ont[,i]+dat$oot[,i],ylim=(c(0,0.5)), pch=19,  cex=0.8, tck = -0.05, bty='l', xlab = " ", ylab=" ")
lines(rept$ent[,i],lty=2)
text(2,0.45,paste(year[dat$it[i]]),cex=1.0)
}
mtext("Trawl length: observed vs predicted", side=3, outer=T)
mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)
plot.new()
plot.new()
plot.new()
plot.new()
plot.new()
# Create observer survey length vs predicted length freq figure 
for (i in 1:dat$nyo){
x <- plot(dat$ono[,i]+dat$ooo[,i],ylim=(c(0,0.6)), pch=19,  cex=0.8, tck = -0.05, bty='l', xlab = " ", ylab=" ")
lines(rept$eno[,i]+rept$eoo[,i],lty=2)
text(2,0.45,paste(year[dat$io[i]]),cex=1.0)
}
mtext("Observer length: observed vs predicted", side=3, line = -23, outer=T)

# Create pot survey length vs predicted length freq figure 
#par(mfrow=c(6,6),mar = c(1.2, 1.2, 1.2, 1.2),oma=c(3,1,2,1))
#par(mar = c(1.2, 1.2, 2, 1.2),oma=c(3,1,2,1))
#par(mfrow=c(6,6),mar = c(1.2, 1.2, 1.2, 1.2),oma=c(3,1,2,1))
#for (i in 1:dat$nyp){
#x <- plot(dat$onp[,i]+dat$oop[,i],ylim=(c(0,0.5)), pch=19,  cex=0.8, tck = -0.05, bty='l', xlab = " ", ylab=" ")
#lines(rept$enp[,i]+rept$eop[,i],lty=2)
#text(2,0.45,paste(year[dat$ip[i]]),cex=1.0)
#}
#mtext("Pot length: observed vs predicted", side=3, outer=T)
mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)

##############################################################################################################
#  Figure for bubble plot
##############################################################################################################
bubbleplot<-function(){
test2$radius <- 1.5*sqrt(abs(test2$fbdata)/pi) 
symbols(test2$years,test2$classes,circles=test2$radius, bg=ifelse(test2$fbdata<0,'white','black'),fg=ifelse(test2$fbdata==0,'white','black'), inches = FALSE, xlim=c(1976,1975+dat$ny-1),xlab = " ", ylab=" ")
}

par(mfrow=c(5,1), mar=c(2,2,2,1),oma=c(3,1,2,1))

# Create commercial length vs predicted length freq figure 
bdata <- dat$onc+dat$ooc - (rept$enc + rept$eoc)
fbdata <- as.vector(bdata)
classes <- rep(1:6,dat$ny)
years <- sort(rep(1976:(1975+dat$ny),6))
test <- data.frame(cbind(years,classes,fbdata))
test2 <- subset(test,((test$years != 1976)&(test$years != 1991)))
bubbleplot()
title(main='Commercial Harvest')

# Create Winter length vs predicted length freq figure 
bdata <- dat$onw+dat$oow - rept$enw - rept$eow
fbdata <- as.vector(bdata)
classes <- rep(1:6,dat$nyw)
years <-1
for (i in 1:dat$nyw){
years <- c(years,rep(year[dat$iw[i]],6))
}
years <- years[-1]
test2 <- data.frame(cbind(years,classes,fbdata))
bubbleplot()
title(main='Winter Pot Survey')

# Create Trawl length vs predicted length freq figure 
bdata <- dat$ont+dat$oot - rept$ent
fbdata <- as.vector(bdata)
classes <- rep(1:6,dat$nyt)
years <-1
for (i in 1:dat$nyt){
years <- c(years,rep(year[dat$it[i]],6))
}
years <- years[-1]
test2 <- data.frame(cbind(years,classes,fbdata))
bubbleplot()
title(main='Trawl Survey')

# Create Observer survey length vs predicted length freq figure 
bdata <- dat$ono+dat$ooo - (rept$eno + rept$eoo) 
fbdata <- as.vector(bdata)
classes <- rep(1:6,dat$nyo)
years <-1
for (i in 1:dat$nyo){
years <- c(years,rep(year[dat$io[i]],6))
}
years <- years[-1]
test2 <- data.frame(cbind(years,classes,fbdata))
bubbleplot()
title(main='Observer Survey')

# Create Pot survey length vs predicted length freq figure 
#bdata <- dat$onp+dat$oop - (rept$enp + rept$eop) 
#fbdata <- as.vector(bdata)
#classes <- rep(1:6,dat$nyp)
#years <-1
#for (i in 1:dat$nyp){
#years <- c(years,rep(year[dat$ip[i]],6))
#}
#years <- years[-1]
#test2 <- data.frame(cbind(years,classes,fbdata))
#bubbleplot()
#title(main='Summer Pot Survey')


mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)


##############################################################################################################
#  Figure for Length Frequency Plot
##############################################################################################################
par(mfrow=c(2,3),mar = c(2, 1.5, 1.5, 1.5),oma=c(4,1,2,1))
lengthclass  <- c(74.0, 84.0, 94.0, 104, 114, 124,134)
# Create commercial length vs predicted length frequencies
cdata <- rowSums((dat$onc+dat$ooc)*dat$tc)
ccdata <- c(cumsum(cdata)/sum(cdata),1)
pdata <- rowSums((rept$enc + rept$eoc)*dat$tc)
cpdata <- c(cumsum(pdata)/sum(pdata),1)
plot(lengthclass, ccdata, type='s', main='Commercial Harvest', xlab = 'Length', ylab='')
lines(lengthclass, cpdata, type='s',lty=5, col='red')

# Create Winter length vs predicted length freq figure 
cdata <- rowSums((dat$onw+dat$oow)*dat$sw)
ccdata <- c(cumsum(cdata)/sum(cdata),1)
pdata <- rowSums((rept$enw + rept$eow)*dat$sw)
cpdata <- c(cumsum(pdata)/sum(pdata),1)
plot(lengthclass, ccdata, type='s', main='Winter Pot Survey', xlab = 'Length', ylab='')
lines(lengthclass, cpdata, type='s', lty=5, col='red')

# Create Trawl length vs predicted length freq figure 
cdata <- rowSums((dat$ont+dat$oot)*dat$tt)
ccdata <- c(cumsum(cdata)/sum(cdata),1)
pdata <- rowSums((rept$ent)*rept$ett)
cpdata <- c(cumsum(pdata)/sum(pdata),1)
plot(lengthclass, ccdata, type='s', main='Trawl Survey', xlab = 'Length', ylab='')
lines(lengthclass, cpdata, type='s', lty=5, col='red')

# Create Observer survey length vs predicted length freq figure 
cdata <- rowSums((dat$ono+dat$ooo)*dat$so)
ccdata <- c(cumsum(cdata)/sum(cdata),1)
pdata <- rowSums((rept$eno + rept$eoo)*dat$so)
cpdata <- c(cumsum(pdata)/sum(pdata),1)
plot(lengthclass, ccdata, type='s', main='Summer Observer', xlab = 'Length', ylab='')
lines(lengthclass, cpdata, type='s', lty=5, col='red')

# Create Pot survey length vs predicted length freq figure 
#cdata <- rowSums((dat$onp+dat$oop)*dat$tp)
#ccdata <- c(cumsum(cdata)/sum(cdata),1)
#pdata <- rowSums((rept$enp + rept$eop)*rept$etp)
#cpdata <- c(cumsum(pdata)/sum(pdata),1)
#plot(lengthclass, ccdata, type='s', main='Summer Pot Survey', xlab = 'Length', ylab='')
#lines(lengthclass, cpdata, type='s', lty=5, col='red')
#mtext("Length frequency: Observed (black) vs. Predicted (red)", side=3,outer=T)
#dev.off()

##############################################################################################################
#  Figure for length stack graph
##############################################################################################################

par(mfrow=c(5,1), mar=c(2,2,2,1),oma=c(3,1,2,1))
color <- c('white', 'gray', 'dimgray', 'darkgray','black','gray90') 
# Create commercial length freq figure 
bdata <- dat$onc+dat$ooc 
barplot(bdata[,1:dat$ny],names.arg = year[1:dat$ny], col =color)
legend('left',legend=c(6,5,4,3,2,1), fill=c('gray90','black','darkgray','dimgray','gray','white') , bg='white')
title(main='Commercial Harvest')

mdata <- matrix(0,6,dat$ny)
# Create Winter length freq figure 
bdata <- dat$onw+dat$oow
for (j in 1:dat$nyw) {
	mdata[,dat$iw[j]] <- bdata[,j]
	}
barplot(mdata,names.arg = year[1:dat$ny], col =color)	
title(main='Winter Pot Survey')

mdata <- matrix(0,6,dat$ny)
# Create Summer Trawl length freq figure 
bdata <- bdata <- dat$ont+dat$oot 
for (j in 1:dat$nyt) {
	mdata[,dat$it[j]] <- bdata[,j]
	}
barplot(mdata,names.arg = year[-dat$ny], col =color)	
title(main='Summter Trawl Survey')

mdata <- matrix(0,6,dat$ny)
# Create Observer freq figure 
bdata <- bdata <- dat$ono+dat$ooo 
for (j in 1:dat$nyo) {
	mdata[,dat$io[j]] <- bdata[,j]
	}
barplot(mdata,names.arg = year[-dat$ny], col =color)	
title(main='Summer Observer Survey')

mdata <- matrix(0,6,dat$ny)
# Create Summer Pot length freq figure 
#bdata <- bdata <- dat$onp+dat$oop 
#for (j in 1:dat$nyp) {
#	mdata[,dat$ip[j]] <- bdata[,j]
#	}
#barplot(mdata,names.arg = year[-dat$ny], col=color)
#legend('left',legend=c(6,5,4,3,2,1), fill=c('gray90','black','darkgray','dimgray','gray','white') , bg='white')
#title(main='Summter Pot Survey')
mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)



dev.off()



















