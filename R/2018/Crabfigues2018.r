############################################################################
#  Figure producing R function rutine
############################################################################
rm(list=ls(all=TRUE))
library(gplots)
library(PBSadmb)
library(mgcv)
admodeldir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2018/Model/')
setwd(admodeldir)

############################################################################
#   1.0 Read source file 
############################################################################
# Read read.R for program use
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/ADMB_rw_data.R')
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Prediction_model_dataset/Crab_model_data.R')

# Crab Julyschedule change routine
#source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/make_July_data.R')

############################################################################
#   1.1 Define input file name and output file name and locations. 
############################################################################
# infun: input file name 
datadir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/')

############################################################################
#   Activate either, based on model you use 
############################################################################
#outfn <- paste(datadir,'ns3n2016R_5mm.dat',sep='')
outfn <- paste(admodeldir,'ns3n2018new4.dat',sep='')
pdffile <- paste(admodeldir,'ns3n2018.pdf',sep='')

############################################################################
#   1.1.1 ADMB program file locations: Choose prefix based on your modeling 
############################################################################
# ADMB model file name
prefix <- paste(admodeldir,'ns3n2018_Feb1_tr_new',sep='')

# fnpar: ADMB parameter ouput file name 
fn <- paste(prefix,'.rep',sep='')
fnpar <- paste(prefix,'.std',sep='')


############################################################################
#   1.2 Red ADMB data into R  
############################################################################
# Read crab data into R List file 
dat<-list()
######## Create ADMB data ##################################################
dat$fyear <- 1976
dat$lyear <- 2017
# The number of length class 
dat$na <- dim(lengthclass)[1]
# The number of recruit length class 
dat$ra <- length(lengthclass[(which(lengthclass$lenclass < 94)),])
# The number of sublegal length class 
dat$sla <- length(lengthclass[(which(lengthclass$lenclass < 104)),])
# The number of trawl survey 
dat$nyt <- length(trawl$Year)
# The number of winter pot survey 
dat$nyw <- length(wntr$Year)
# The number of observer survey 
dat$nyo <- length(ob.yrs)
# Row number of tag data 
dat$ntag1 <- dim(t1w)[1]
# Row number of tag data 
dat$ntag2 <- dim(t2w)[1]
# summer commercial fishery event indices: 1)year large trawl fishery ended, 2) escape mechanism installed, 3) year summer commrecial fishery CW>5 inch crab was accepted by buyers
dat$scy <- c(1992,1993,2005)
# year NOAA survey ended 
dat$qyear <- 1992
# mid length of the smallest length class
dat$slm <- min.s+(inc.s-1)/2  
# interval of each length group
dat$slt <- inc.s  
# Natural mortality
dat$M <- 0.18
# Natural mortality multiplier for the last length class
dat$msn <- c(length(lengthclass[(which(lengthclass$lenclass >= 124)),]),length(lengthclass[(which(lengthclass$lenclass < 84)),]) )
# Maximum effective sample size for (1) trawl and (2) commercial and observer
dat$maxs <- c(20,10)
# proporion of sample size for (1) trawl survey  and (2) winter pot, summer com, observer
dat$efn <- c(0.5,0.1)
# Discards handling mortality for (1) summer and (2) winter
dat$hm <- c(0.2,0.2)
# The proportion of legals by length group
dat$lg <- p.legal
# Mean weight by length class
dat$wm <- lwp
# Years trawl survey conducted 
dat$it <- trawl$Year
# Trawl survey aubndace 
dat$tt <- trawl$Abundance
# Proportion of Com catch occred before mid survey date 
dat$pct <- trawl$pct
# Mid point of trawl survey from July 01 
dat$yt <- trawl$yt
# Trawl survey CV
dat$cv <- trawl$CV
# Trawl survey newshell length 
dat$nont <- as.matrix(trwl.size[trwl.size$Shell==1,-c(1:2)])
# Trawl survey oldshell length 
dat$noot <- as.matrix(trwl.size[trwl.size$Shell==2,-c(1:2)])
# Winter Pot survey years conducted 
dat$iw <- wntr$Year
# Winter Pot survey CPUE
dat$wcpue <- wntr$CPUE
# Winter Pot survey newshell length
dat$nonw <- as.matrix(wntr.size[wntr.size$Shell==1,-c(1:2)])
# Winter Pot survey oldshell length
dat$noow <- as.matrix(wntr.size[wntr.size$Shell==2,-c(1:2)])
# Summer Commercial Catch 
dat$tc <- Harvest$S.Com
# Summer Commercial Pot lift effort  
dat$te <- Harvest$Pot.Lift
# Summer Commercial Standardized CPUE  
dat$stcpue <- Harvest$ST.CPUE
# Summer Commercial Standardized Starndard error  
dat$secpue <- Harvest$ST.SE
# Mid point of Commercial harvest from July 01 
dat$ys <- Harvest$ys
# Summer Commercial newshell length
dat$nonc <- as.matrix(com.size[com.size$Shell==1,-c(1:2)])
# Summer Commercial oldshell length
dat$nooc <- as.matrix(com.size[com.size$Shell==2,-c(1:2)])
# Winter Commercial fishery
dat$twc <- Harvest$W.Com
# Winter Subsistence fishery Retained
dat$tws <- Harvest$W.Sub.R
# Winter Subsistence fishery Caught Total 
dat$twst <- Harvest$W.Sub.T
# Observer survey years conducted 
dat$io <- ob.yrs
# Summer Commercial Observer discards newshell length
dat$nono <- as.matrix(obs.size[obs.size$Shell==1,-c(1:2)])
# Summer Commercial Observer discards oldshell length
dat$nooo <- as.matrix(obs.size[obs.size$Shell==2,-c(1:2)])
# Tag recovery data first period
dat$tag_recov1 <- as.matrix(t1w)
# Tag recovery data second period
dat$tag_recov2 <- as.matrix(t2w)

############################################################################
#   1.3 Update data  
############################################################################
########  Add weights to Summer pot survey likelihood ######################
dat$SDRec <- 0.5
dat$SDW <- 0.3
dat$M_phase <- -1
dat$rmol_phase <- -1
dat$lamc <- 1
dat$lawp <- 1
dat$latag <- 0.5
dat$nst <- 1
dat$nsc <- 1
dat$ssc <- -1
dat$pwh <- 0.16


############################################################################	

############################################################################	
		
makeADdata(dat,outfn)
outfn <- paste(admodeldir,'ns3n2018new4nost.dat',sep='')

# Run the ADMB 
argvec <- paste(' -ind ',outfn,sep="")
runAD(prefix=prefix, argvec = argvec, logfile=TRUE, add = FALSE, verbose=TRUE)

# Run ADMB with MCMC 
#argvec <- paste(' -ind ',outfn,sep="", ' -mcmc 1000000 -mcsave 1000')
#runAD(prefix=prefix, argvec = argvec, logfile=TRUE, add = FALSE, verbose=TRUE)
#argvec <- paste(' -ind ',outfn,sep="", ' -mceval')
#runAD(prefix=prefix, argvec = argvec, logfile=TRUE, add = FALSE, verbose=TRUE)

# Read output file
rept <- reptoRlist(fn)
#rept
par1 <- read.table(fnpar, sep="",header=T)
pars <- data.frame(t(par1[,3]))
names(pars) <- t(par1[,2])
rept$f
rept$tf
par1
tr <- rept$mgr
b <- c(1,0,0,0,0,0,0,0)
s <- c(68.5,78.5,88.5,98.5,108.5,118.5,128.5,138.5)
m <- matrix(0,8,20)
for(i in 1:20) {
if(i==1){b1 <- b}
m[,i] <- t(tr)%*%b1
b1 <- m[,i] 
}
m1 <- cbind(b,m)

for(i in 1:8){
if(i ==1){plot(m1[i,],type='l')}
else{
lines(m1[i,])
}
}
plot(colSums(m1*s),type='l')



############################################################################
#  Calculate Sample numbers 
############################################################################
# Calculate total 
  st <- colSums(dat$nont)+colSums(dat$noot)
  sw <- colSums(dat$nonw)+colSums(dat$noow)
  sc <- colSums(dat$nonc)+colSums(dat$nooc)
  so <- colSums(dat$nono)+colSums(dat$nooo)

  ont <- dat$nont/st[col(dat$nont)]
  oot <- dat$noot/st[col(dat$noot)]
  onw <- dat$nonw/sw[col(dat$nonw)]
  oow <- dat$noow/sw[col(dat$noow)]
  onc <- dat$nonc/sc[col(dat$nonc)]
  ooc <- dat$nooc/sc[col(dat$nooc)]
  ono <- dat$nono/so[col(dat$nono)]
  ooo <- dat$nooo/so[col(dat$nooo)]

############################################################################
#  Calculate Sample numbers 
############################################################################
year <- seq(dat$fyear,dat$lyear)
year.1 <- c(year,dat$lyear+1)
# Define number of years 
ny <- length(year)
length.c  <- seq(min.s,max.s+inc.s,inc.s)
sizeclass <- length.c + (inc.s-1)/2
sizeclass[length(sizeclass)] <- sizeclass[length(sizeclass)-1]


############################################################################
#  Figures 
############################################################################
windows(width=20,height=20,record=TRUE)
 
############################################################################
#  Plot Annual Molting probability 
############################################################################
if(dat$rmol_phase==1)
{
par(mfrow=c(6,7),mar = c(2, 1.5, 1.5, 1.5),oma=c(4,1,2,1))
for(i in 1:ny){
moltp <- moltp <- c(rept$molp[,i],rept$molp[dat$na,i])
plot(length.c, moltp,type='s', ylim=c(0,1),,main = year[i])
}
mtext("Annual Molting probability", side=3, outer=T)
}
############################################################################
#  Plot Selectivity
############################################################################
par(mfrow=c(2,3),mar = c(2, 1.5, 1.5, 1.5),oma=c(4,1,2,1))
# Molting Probability 
moltp <- c(rept$molp[,1],rept$molp[dat$na,1])
if(dat$rmol_phase==-1) plot(length.c, moltp,type='s',ylim=c(0,1), main='Molting Probability')
# Trawl Selectivity
# NOAA or 1 trawl selectivity
trawlpn <- c(rept$selt1[1,],rept$selt1[1,dat$na])
plot(length.c, trawlpn,type='s',ylim=c(0,1))
if(dat$nst==2) title(main='NOAA Trawl Selectivity') else title(main='Trawl Selectivity')
# ADFG or 2 trawl selectivity
trawlpa <- c(rept$selt1[2,],rept$selt1[2,dat$na])
if(dat$nst==2) plot(length.c, trawlpa,type='s', main='ADFG Trawl Selectivity', ylim=c(0,1))
# Winter Pot Selectivity
wpotp <- c(rept$selw,rept$selw[dat$na])
plot(length.c, wpotp,type='s', main='Winter pot Selectivity', ylim=c(0,1))
# Summer Commercial Selectivity
sc1p <- c(rept$selc1[1,],rept$selc1[1,dat$na])
plot(length.c, sc1p,type='s',ylim=c(0,1))
if(dat$nsc==2) title(main='Commercial Selectivity 76-92') else title(main='Commercial Selectivity')
sc2p <- c(rept$selc1[2,],rept$selc1[2,dat$na])
if(dat$nsc==2) plot(length.c, sc2p,type='s', main='Commercial Selectivity 93-15', ylim=c(0,1))
#sc3p <- selc(pars$log_sc5)
#plot(lengthclass, sc3p,type='s', main='Commercial 93-14 Selectivity', ylim=c(0,1))


############################################################################
#  Plot residual QQ plot for survey
############################################################################
############################################################################
#  Calculate RMSE and Residual plots
############################################################################
trawl1 <- ifelse(is.null(pars$qtno),1,pars$qtno)*(rept$ett[1:2])
trawl2 <- ifelse(is.null(pars$qtad),1,pars$qtad)*(rept$ett[3:dat$nyt])
etrawl <- c(trawl1,trawl2)
RMSE.trawl <- sqrt(sum((log(dat$tt)-log(etrawl))^2)/length(dat$tt))
#RMSE.pot <- sqrt(sum((log(dat$tp)-log(rept$etp1+rept$etp2))^2)/length(dat$tp))
RMSE.cpue <- sqrt(sum((log(dat$stcpue[-c(1,13,16,ny)])-log(rept$ecpue[-c(1,13,16,ny)]))^2)/length(dat$stcpue[-c(1,13, 16,ny)]))
RMSE.trawl
#RMSE.pot
RMSE.cpue
resd.trawl <- dat$tt - etrawl
resd.pot <- dat$tp - (rept$etp1+rept$etp2)
resd.cpue <- dat$stcpue[-c(1,13,16,ny)] - rept$ecpue[-c(1,13,16,ny)]

par(mfrow=c(3,3),mar = c(2, 2, 2, 2),oma=c(3,1,2,1))
hist(resd.trawl, main='Trawl survey')
qqnorm(resd.trawl)
qqline(resd.trawl)
plot((etrawl),resd.trawl,main='predicted vs. residual')
legend('topright',paste('RMSE',round(RMSE.trawl,3)),bty='n',xjust=0)
abline(h=0)
hist(resd.cpue,, main='Commercial CPUE')
qqnorm(resd.cpue)
qqline(resd.cpue)
plot(rept$ecpue[-c(1,13,16,ny)],resd.cpue,main='predicted vs. residual')
legend('topright',paste('RMSE',round(RMSE.cpue,3)),bty='n',xjust=0)
abline(h=0)
#hist(resd.pot,main='Summer pot survey Residuals')
#qqnorm(resd.pot)
#qqline(resd.pot)
#plot((rept$etp1+rept$etp2),resd.pot)
#abline(h=0)
mtext("Residuals Histogram, Q-Q Plot, Predicted vs. Residual", side=3, outer=T)

############################################################################
#  Calculate and Plot Effective sample size 
############################################################################
#trawl survey effective sample size
ot <- ont+oot
tefn <-colSums(rept$ent*(1-rept$ent))/colSums((ot-rept$ent)^2)
tefnn <- ifelse(st*dat$efn[1]>dat$maxs[1],dat$maxs[1],st*dat$efn[1])
 
#Summer pot survey effective sample size
#spp <- rept$enp+rept$eop
#sp <- dat$onp+dat$oop
#spefn <-colSums(spp*(1-spp))/colSums((spp-sp)^2)
#spefnn <- ifelse(dat$sp*0.5>dat$maxss,dat$maxss,dat$sp*0.5)

#commercial fishery effective sample size
scp <- rept$enc+rept$eoc
os <- onc+ooc
scefn <-colSums(scp*(1-scp))/colSums((scp-os)^2)
scefn <-scefn[-c(1,16)]
scefnn <- ifelse(sc*dat$efn[2]>dat$maxs[2],dat$maxs[2],sc*dat$efn[2])[-c(1,16)]

#Winter Pot survey effective sample size
wpp <- rept$enw+rept$eow
wp <- onw+oow
wpefn <-colSums(wpp*(1-wpp))/colSums((wpp-wp)^2)
wpefnn <- ifelse(sw*dat$efn[2]>dat$maxs[2],dat$maxs[2],sw*dat$efn[2])

#Tag recovery survey effective sample size;
# Create empty matrix
tagrecap11 <- matrix(0,dat$na,dat$na)
tagrecap21 <- matrix(0,dat$na,dat$na)
tagrecap31 <- matrix(0,dat$na,dat$na)
tagrecap12 <- matrix(0,dat$na,dat$na)
tagrecap22 <- matrix(0,dat$na,dat$na)
tagrecap32 <- matrix(0,dat$na,dat$na)
tagrecap1 <- matrix(0,dat$na,dat$na)
tagrecap2 <- matrix(0,dat$na,dat$na)
tagrecap3 <- matrix(0,dat$na,dat$na)    
# fill tag recovery data 
for(i in 1:dat$ntag1) 
  {
        tagrecap11[dat$tag_recov1[i,1],dat$tag_recov1[i,2]] <- dat$tag_recov1[i,3]
        tagrecap21[dat$tag_recov1[i,1],dat$tag_recov1[i,2]] <- dat$tag_recov1[i,4]
        tagrecap31[dat$tag_recov1[i,1],dat$tag_recov1[i,2]] <- dat$tag_recov1[i,5]
  }
# fill tag recovery data 
for(i in 1:dat$ntag2) 
  {
        tagrecap12[dat$tag_recov2[i,1],dat$tag_recov2[i,2]] <- dat$tag_recov2[i,3]
        tagrecap22[dat$tag_recov2[i,1],dat$tag_recov2[i,2]] <- dat$tag_recov2[i,4]     
        tagrecap32[dat$tag_recov2[i,1],dat$tag_recov2[i,2]] <- dat$tag_recov2[i,5]  
  }
# Combine two matrices   
tagrecap1 <- tagrecap11+tagrecap12
tagrecap2 <- tagrecap21+tagrecap22
tagrecap3 <- tagrecap31+tagrecap32

  
# calculate input sample size (sum of each row frequency)
nrecap <- c(rowSums(tagrecap1),rowSums(tagrecap2), rowSums(tagrecap3),rowSums(tagrecap12), rowSums(tagrecap22),rowSums(tagrecap32))
#nrecap <- c(rowSums(tagrecap1),rowSums(tagrecap2), rowSums(tagrecap3))

# change tag recovery to proportion
tagrecap11 <- tagrecap11/rowSums(tagrecap11)
tagrecap21 <- tagrecap21/rowSums(tagrecap21)
tagrecap31 <- tagrecap31/rowSums(tagrecap31)
tagrecap12 <- tagrecap12/rowSums(tagrecap12)
tagrecap22 <- tagrecap22/rowSums(tagrecap22)
tagrecap32 <- tagrecap32/rowSums(tagrecap32)

tagrecap1 <- tagrecap1/rowSums(tagrecap1)
tagrecap2 <- tagrecap2/rowSums(tagrecap2)
tagrecap3 <- tagrecap3/rowSums(tagrecap3)

# replace NaN to 0 
tagrecap1[is.nan(tagrecap1)] <- 0
tagrecap2[is.nan(tagrecap2)] <- 0
tagrecap3[is.nan(tagrecap3)] <- 0
tagrecap11[is.nan(tagrecap11)] <- 0
tagrecap21[is.nan(tagrecap21)] <- 0
tagrecap31[is.nan(tagrecap31)] <- 0
tagrecap12[is.nan(tagrecap12)] <- 0
tagrecap22[is.nan(tagrecap22)] <- 0
tagrecap32[is.nan(tagrecap32)] <- 0

# calculate effective sample size 
esn1 <- rowSums(rept$egr1*(1-rept$egr1))/rowSums((rept$egr1-tagrecap1)^2)
esn2 <- rowSums(rept$egr2*(1-rept$egr2))/rowSums((rept$egr2-tagrecap2)^2)
esn3 <- rowSums(rept$egr3*(1-rept$egr3))/rowSums((rept$egr3-tagrecap3)^2)
esn12 <- rowSums(rept$egr12*(1-rept$egr12))/rowSums((rept$egr12-tagrecap12)^2)
esn22 <- rowSums(rept$egr22*(1-rept$egr22))/rowSums((rept$egr22-tagrecap22)^2)
esn32 <- rowSums(rept$egr32*(1-rept$egr32))/rowSums((rept$egr32-tagrecap32)^2)
enrecap <- c(esn1, esn2, esn3, esn12, esn22, esn32)
#enrecap <- c(esn1, esn2, esn3)

# replace NaN to 0 
enrecap[is.nan(enrecap)] <- 0


# Observer survey effective sample size
opp <- rept$eno+rept$eoo
op <- ono+ooo
opefn <-colSums(opp*(1-opp))/colSums((opp-op)^2)
opefnn <- ifelse(so*dat$efn[2]>dat$maxs[2],dat$maxs[2],dat$so*dat$efn[2])

############################################################################
#  Plot Effective sample size 
############################################################################
par(mfrow=c(5,3),mar = c(3,3,2, 2),oma=c(3,1,2,1))

# Trawl Survey Figure
tmain <- 'Trawl survey'
# plot Implied Frequency
hist(tefn,breaks=10, main=tmain)
title(xlab = 'Implied', line = 2)
title(ylab = 'Frequency', line = 2)
abline(v=mean(tefn),lwd=2)
# plot Input vs. Implied
plot(tefnn,tefn, main=tmain, xlim = c(0,dat$maxs[1]),ylim=c(0,max(tefn)))
abline(lm(tefn~ tefnn -1),lty=2)
abline(0,1)
title(xlab = 'Input', line = 2)
title(ylab = 'Implied', line = 2)
# Plot Implied vs. year
plot(dat$it, tefn,main=tmain)
title(xlab = 'Year', line = 2)
title(ylab = 'Implied', line = 2)

# Commercial Catch Figure
tmain <- 'Commercial Catch'
# plot Implied Frequency
hist(scefn,breaks=10, main=tmain)
abline(v=mean(scefn),lwd=2)
title(xlab = 'Implied', line = 2)
title(ylab = 'Frequency', line = 2)
# plot Input vs. Implied
plot(scefnn, scefn, main=tmain, xlim = c(0,dat$maxs[2]),ylim=c(0,max(scefn)))
lines(scefn~scefnn-1,lty=2)
abline(0,1)
title(xlab = 'Input', line = 2)
title(ylab = 'Implied', line = 2)
# Plot Implied vs. year
plot(year[-c(1,16)], scefn,main=tmain)
title(xlab = 'Year', line = 2)
title(ylab = 'Implied', line = 2)


# Winter Pot Survey Figure
tmain <- 'Winter pot survey'
# Plot Implied Frequency
hist(wpefn,breaks=5, main=tmain)
abline(v=mean(wpefn),lwd=2)
title(xlab = 'Implied', line = 2)
title(ylab = 'Frequency', line = 2)
# plot Input vs. Implied
plot(wpefnn, wpefn, main=tmain, xlim = c(0,dat$maxs[2]),ylim=c(0,max(wpefn)))
abline(lm(wpefn~wpefnn-1),lty=2)
abline(0,1)
title(xlab = 'Input', line = 2)
title(ylab = 'Implied', line = 2)
# Plot Implied vs. year
plot(dat$iw, wpefn, main=tmain)
title(xlab = 'Year', line = 2)
title(ylab = 'Implied', line = 2)

# Observer Survey Figure
tmain <- 'Observer survey'
# Plot Implied Frequency
hist(opefn,breaks=10, main=tmain)
abline(v=mean(opefn),lwd=2)
title(xlab = 'Implied', line = 2)
title(ylab = 'Frequency', line = 2)
# plot Input vs. Implied
plot(opefnn, opefn, main=tmain, xlim = c(0,dat$maxs[2]),ylim=c(0,max(opefn)))
abline(lm(opefn~opefnn-1),lty=2)
abline(0,1)
title(xlab = 'Input', line = 2)
title(ylab = 'Implied', line = 2)
# Plot Implied vs. year
plot(dat$io, opefn, main=tmain)
title(xlab = 'Year', line = 2)
title(ylab = 'Implied', line = 2)

# Tag recovery Figure
tmain <- 'Tag recovery'
# Plot Implied Frequency
hist(enrecap,breaks=10, main=tmain)
abline(v=mean(enrecap),lwd=2)
title(xlab = 'Implied', line = 2)
title(ylab = 'Frequency', line = 2)
# plot Input vs. Implied
plot(nrecap, enrecap, main=tmain)
abline(lm(enrecap~nrecap-1),lty=2)
abline(0,1)
title(xlab = 'Input', line = 2)
title(ylab = 'Implied', line = 2)

# Summer Pot Survey Figure
#hist(spefn,breaks=10, main='Summer pot survey')
#abline(v=c(dat$maxss),lwd=2)
#plot(spefnn,spefn, main='Summer pot survey', xlim = c(0,dat$maxss),ylim=c(0,max(spefn)))
#abline(lm(spefn ~ spefnn -1),lty=2)
#lines(c(0,dat$maxss),c(0,dat$maxss))
#plot(dat$ip+1975, spefn)
mtext("Effective sample size", side=3, outer=T)


############################################################################
#  Figure for trends.
############################################################################
par(mfrow=c(1,1), mar=c(5,5,4,4))
# Figure 1. Total estimated abundance 
# Plot trawl survey data with CI
qtno <- c(rep(pars$qtno,6),rep(1,dat$nyt-6))
t-t <- dat$tt/qtno
\ciw <- 2*tt*dat$cv/1000-*
ciwl <- exp(log(tt/1000)*(1+2*dat$cv))
ciwu <- exp(log(tt/1000)*(1-2*dat$cv))
plot(dat$it,tt/1000, ylab = 'Crab Abundance (million)', xlab = 'Year', ylim= c(0,8),bty='l')
arrows(dat$it,y0=ciwu,y1=ciwl,code=3,length=0.025, angle =90)
points(dat$it, tt/1000, pch=21, bg='white')
#points(year[dat$ip], dat$tp/1000, pch= 19)
# Model estiamted trawl survey abundance
points(dat$it,etrawl/1000,pch=19, col=4)
#lines(year,total.abundance[1:ny],lwd=2)
legend("topright", pch = c(1,20,20), col= c(1,4,1), legend = c('Observed','Predicted'))
title(main='Trawl survey crab abundance')


# Figure 2. Total, Recruites, Legals estimated abundance 
#: Add projected number 
#rept$nps[1:5,ny+1] <- rept$npw[1:5,ny+1]*exp(-0.417*dat$M)
#rept$nps[6,ny+1] <- rept$npw[6,ny+1]*exp(-0.417*dat$M*dat$ms)
#rept$ops[1:5,ny+1] <- rept$opw[1:5,ny+1]*exp(-0.417*dat$M)
#rept$ops[6,ny+1] <- rept$opw[6,ny+1]*exp(-0.417*dat$M*dat$ms)

total.abundance <- colSums(rept$npw+rept$opw)/1000
plot(year.1, total.abundance, ylab = 'Abundance (million crabs)', typ = 'l', lwd=2, xlab = 'Year',ylim= c(0,10),bty='l')
legal.abundance <- colSums(rept$npw[(dat$sla+1):dat$na,]+rept$opw[(dat$sla+1):dat$na,])
lines(year.1, legal.abundance/1000, lty=2)
recruit.abundance <- colSums(rept$npw[1:dat$ra,]+rept$opw[1:dat$ra,])
lines(year.1, recruit.abundance/1000, lty=3)
legend("topright",lty = c(1,2,3), legend = c('total','legal','recruits'))
title(main='Modeled crab abundance Feb 01')


# Figure 5. MMB data 
mmb <- rept$MMB[-length(rept$MMB)]
plot(year, mmb/1000, ylab = 'MMB (million lb)', lwd=2, typ = 'l', xlab = 'Year', bty='l')
# Bmsy is from 1980
bmsy <- mean(rept$MMB[-c(1:4)])
abline(h=bmsy/1000, lwd=2, lty=2)
points(dat$lyear+1,rept$MMB[length(rept$MMB)]/1000,pch=19,col=1)
# Calculate ratio of MMB to Bmsy
# p.y is projected column length
p.y <- length(rept$MMB)
cofl <- rept$MMB[length(rept$MMB)]/bmsy
# calculate Bmsy adjusted coefficient
fl <- ifelse(cofl>1,1,ifelse(0.25<cofl&cofl<=1,(cofl-0.1)/0.9,0))
# calculate FOF:
FOFL <- rept$M*fl
# Calculate Legal biomass = projected winter*fishery selectivity
#                           *prop legal*weight 
Legal.b <- with(rept,npw[,p.y]+opw[,p.y])
			*dat$lg*rept$selc[,p.y-1]*dat$wm
Legal.s.b <- with(rept,npw[,p.y]+opw[,p.y])*(1-dat$lg)*rept$selc[,p.y-1]*dat$wm			
# Calculate adjusted FOFL 
par1 <- exp(-(FOFL +0.42*rept$M))
par2 <- 1-exp(-0.42*rept$M)
pw <- dat$pwh
FOFL.a <- 1 - par1 - par2*(1-pw*par1)/(1-pw*par2)
OFL.r <- sum(Legal.b*FOFL.a)
OFL.nr <- sum(Legal.s.b*(FOFL.a)*dat$hm[1])



#OFLrb <- (1-exp(-FOFL))*pars$last_legal/1000
title(main='MMB Feb 01')
#tex <- c(paste('BMSY ',round(bmsy/1000,3),' mil.lb'),paste('MMB ',round(pars$last_mmb/1000,3),' mil.lb'),paste('Legal B ',round(pars$last_legal/1000,3),' mil.lb'),paste('FOFL ',round(FOFL,3)), paste('OFL ',round(OFLrb,3),' mil.lb'),paste('ABC ',round(ABC,3),' mil.lb'))
tex <- c(paste('BMSY ',round(bmsy/1000,3),' mil.lb'),paste('MMB ',round(rept$MMB[length(rept$MMB)]/1000,3),' mil.lb'),paste('Legal B ',round(rept$Legalb[length(rept$Legalb)]/1000,3),' mil.lb'))
legend("topright",tex,bty='n',xjust=0) 


# Figure 3. cpue data 
dat$stcpue[dat$stcpue==0]<-NA
ciw <- 2*dat$secpue
ciw2 <- 2*dat$stcpue*sqrt(exp(log((dat$secpue/dat$stcpue)^2+1)+pars$advar)-1)
ciw[ciw==0]<-NA
plot(year, dat$stcpue, ylab = 'CPUE',  xlab = 'Year', ylim=c(0,7), bty='l')
arrows(year,y0=dat$stcpue+ciw2,y1=dat$stcpue-ciw2,code=0,col='red')
arrows(year,y0=dat$stcpue+ciw,y1=dat$stcpue-ciw,code=3,length=0.025, angle =90)
points(year, dat$stcpue, pch=21, bg='white')
rept$ecpue[rept$ecpue==0]<-NA
lines(year,rept$ecpue,lwd=2)
legend("topright",lty =c(1,0), lwd = c(2,0), pch=c(NA,19), legend = c('Predicted','Observed'))
title(main='Summer commercial standardized cpue')
 

# Figure 4. Catch data 
par(mar=c(5,4,4,5))
tcatch <- dat$tc+dat$twc+dat$tws
hrate <- tcatch/rept$Legal[-(ny+1)]
plot(year, tcatch/1000, ylab = 'Total Catch (million)', lwd=2, typ = 'l', xlab = 'Year', bty='u')
par(new=TRUE)
plot(year,hrate,lty=2, lwd=2, typ = 'l',bty='l', xaxt='n', yaxt='n', xlab='',ylab='', ylim = c(0,max(hrate)))
axis(4)
mtext("Estimated harvest rate",side=4,line=3)
legend("topright", lty = c(1,2), lwd = 2, legend=c('Total Catch', 'Estimated Harvest Rate'),bty='n',xjust=0)
title(main='Total catch & Harvest rate')


############################################################################
#  Figure for size proportion 
############################################################################
len.c <- seq(min.s,max.s,inc.s) 
# Create commercial length vs predicted length freq figure 
par(mfrow=c(6,7),mar = c(1, 1, 1, 1),oma=c(3,2,2,1),mgp=c(3,.3,0))
for (i in 2:15){
plot(len.c,onc[,i]+ooc[,i],ylim=c(0,0.7),pch=19, cex=0.8,tck = -0.05, bty='l')
lines(len.c,rept$enc[,i]+rept$eoc[,i],lty=2)
text(len.c[3],0.6,year[i],cex=1.0)
}
for (i in 17:length(year)){
plot(len.c,onc[,i]+ooc[,i],ylim=c(0,0.7),pch=19, cex=0.8,tck = -0.05, bty='l')
lines(len.c,rept$enc[,i]+rept$eoc[,i],lty=2)
text(len.c[3],0.6,year[i],cex=1.0)
}
mtext("commercial harvest length: observed vs predicted", side=3, outer=T)
#mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)

# Create winter length vs predicted length freq figure 
par(mfrow=c(6,7),mar = c(1, 1, 1, 1),oma=c(3,1,2,1))
for (i in 1:dat$nyw){
plot(len.c,onw[,i]+oow[,i],ylim=c(0,0.5), pch=19, cex=0.8,tck = -0.05,bty='l', xlab = " ", ylab=" ")
lines(len.c,rept$enw[,i]+rept$eow[,i],lty=2)
text(len.c[3],0.45,dat$iw[i],cex=1.0)
}
mtext("Winter pot length: observed vs predicted", side=3, outer=T)
#mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)


# Create trawl length vs predicted length freq figure 
par(mfrow=c(6,7),mar = c(1, 1, 1, 1),oma=c(3,1,2,1))
for (i in 1:dat$nyt){
x <- plot(len.c,ont[,i]+oot[,i],ylim=(c(0,0.5)), pch=19,  cex=0.8,tck = -0.05, bty='l', xlab = " ", ylab=" ")
lines(len.c,rept$ent[,i],lty=2)
text(len.c[3],0.45,dat$it[i],cex=1.0)
}
mtext("Trawl length: observed vs predicted", side=3, outer=T)
#mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)
plot.new()
plot.new()
plot.new()
plot.new()
plot.new()
plot.new()
# Create observer survey length vs predicted length freq figure 
for (i in 1:dat$nyo){
x <- plot(len.c,ono[,i]+ooo[,i],ylim=(c(0,0.6)), pch=19, cex=0.8,tck = -0.05, bty='l', xlab = " ", ylab=" ")
lines(len.c,rept$eno[,i]+rept$eoo[,i],lty=2)
text(len.c[3],0.55,dat$io[i],cex=1.0)
}
mtext("Discards length: observed vs predicted", side=3, line = -30, outer=T)

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
#mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)

############################################################################
#  Tag data 
############################################################################
       
par(mfrow=c(6,(dat$na-1)/1),mar = c(2, 1.2, 1.2, 1.2),oma=c(3,1,2,1),mgp=c(0,0.25,0))
for (i in 1:(dat$na-1)){
x <- plot(len.c,tagrecap1[i,],ylim=c(0,1),pch=19, cex=0.8, tck = -0.05, bty='l',ylab = '', xlab = '')
lines(len.c,rept$egr1[i,],lty=2)
text(len.c[4],0.9,paste(len.c[i],'mm'),cex=1.0)
}
mtext("Tag recovery data observed vs predicted", cex = 0.9,side=3, outer=T)
mtext("Recovery after 1 year", cex = 0.9, side=3, line = -1.5, outer=T)
for (i in 1:(dat$na-1)){
x <- plot(len.c,tagrecap2[i,],ylim=c(0,1),pch=19, cex=0.8, tck = -0.05, bty='l',ylab = '', xlab = '')
lines(len.c,rept$egr2[i,],lty=2)
text(len.c[4],0.9,paste(len.c[i],'mm'),cex=1.0)
}
mtext("Recovery after 2 years", cex = 0.9, side=3, line = -10.5,  outer=T)
for (i in 1:(dat$na-1)){
x <- plot(len.c,tagrecap3[i,],ylim=c(0,1),pch=19, cex=0.8, tck = -0.05, bty='l',ylab = '', xlab = '')
lines(len.c,rept$egr3[i,],lty=2)
text(len.c[4],0.9,paste(len.c[i],'mm'),cex=1.0)
}
mtext("Recovery after 3 years", cex = 0.9,  side=3, line = -20.5, outer=T)



#for (i in 1:(dat$na-1)){
#x <- plot(len.c,tagrecap12[i,],ylim=c(0,1),pch=19, cex=0.8, tck = -0.05, bty='l',ylab = '', xlab = '')
#lines(len.c,rept$egr12[i,],lty=2)
#text(len.c[3],0.9,paste(len.c[i],'mm'),cex=1.0)
#}
#mtext("2008-2015: Recovery after 1 year", cex = 0.9,  side=3, line = -30.5, outer=T)
#for (i in 1:(dat$na-1)){
#x <- plot(len.c,tagrecap22[i,],ylim=c(0,1),pch=19, cex=0.8, tck = -0.05, bty='l',ylab = '', xlab = '')
#lines(len.c,rept$egr22[i,],lty=2)
#text(len.c[3],0.9,paste(len.c[i],'mm'),cex=1.0)
#}
#mtext("2008-2015: Recovery after 2 years", cex = 0.9,  side=3, line = -40.5, outer=T)
#for (i in 1:(dat$na-1)){
#x <- plot(len.c,tagrecap32[i,],ylim=c(0,1),pch=19, cex=0.8, tck = -0.05, bty='l',ylab = '', xlab = '')
#lines(len.c,rept$egr32[i,],lty=2)
#text(len.c[3],0.9,paste(len.c[i],'mm'),cex=1.0)
#}
#mtext("2008-2015: Recovery after 3 years", cex = 0.9,  side=3, line = -50.5, outer=T)
#mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1, cex = 0.9,outer=T)

############################################################################
#  Figure for bubble plot
############################################################################
bubbleplot<-function(){
test2$radius <- 1.5*sqrt(abs(test2$fbdata)/pi) 
symbols(test2$years,test2$classes,circles=test2$radius, bg=ifelse(test2$fbdata<0,'white','black'),fg=ifelse(test2$fbdata==0,'white','black'), inches = FALSE, xlim=c(dat$fyear,dat$lyear),xlab = " ", ylab=" ")
}

par(mfrow=c(5,1), mar=c(2,2,2,1),oma=c(3,1,2,1))

# Create commercial length vs predicted length freq figure 
bdata <- onc+ooc - (rept$enc + rept$eoc)
fbdata <- as.vector(bdata)
classes <- rep(len.c,length(year))
years <- sort(rep(dat$fyear:dat$lyear,dat$na))
test <- data.frame(cbind(years,classes,fbdata))
test2 <- subset(test,((test$years != 1976)&(test$years != 1991)))
bubbleplot()
title(main='Commercial Harvest')

# Create Winter length vs predicted length freq figure 
bdata <- onw+oow - rept$enw - rept$eow
fbdata <- as.vector(bdata)
classes <- rep(len.c,dat$nyw)
years <-1
for (i in 1:dat$nyw){
years <- c(years,rep(dat$iw[i],dat$na))
}
years <- years[-1]
test2 <- data.frame(cbind(years,classes,fbdata))
bubbleplot()
title(main='Winter Pot Survey')

# Create Trawl length vs predicted length freq figure 
bdata <- ont+oot - rept$ent
fbdata <- as.vector(bdata)
classes <- rep(len.c,dat$nyt)
years <-1
for (i in 1:dat$nyt){
years <- c(years,rep(dat$it[i],dat$na))
}
years <- years[-1]
test2 <- data.frame(cbind(years,classes,fbdata))
bubbleplot()
title(main='Trawl Survey')

# Create Observer survey length vs predicted length freq figure 
bdata <- ono+ooo - (rept$eno + rept$eoo) 
fbdata <- as.vector(bdata)
classes <- rep(len.c,dat$nyo)
years <-1
for (i in 1:dat$nyo){
years <- c(years,rep(dat$io[i],dat$na))
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

#mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)

############################################################################
#  Figure for length stack graph
############################################################################

par(mfrow=c(5,1), mar=c(2,2,2,1),oma=c(3,1,2,1))
color <- c('white', 'gray', 'dimgray', 'darkgray','black','gray90') 
# Create commercial length freq figure 
bdata <- onc+ooc 
barplot(bdata,names.arg = year, col =gray.colors(dat$na))
title(main='Commercial Harvest')

mdata <- matrix(0,dat$na,ny)
# Create Winter length freq figure 
bdata <- onw+oow
for (j in 1:dat$nyw) {
	mdata[,dat$iw[j]-dat$fyear+1] <- bdata[,j]
	}
barplot(mdata,names.arg = year, col =gray.colors(dat$na))
title(main='Winter Pot Survey')

mdata <- matrix(0,dat$na,ny)
# Create Summer Trawl length freq figure 
bdata <- bdata <- ont+oot 
for (j in 1:dat$nyt) {
	mdata[,dat$it[j]-dat$fyear+1] <- bdata[,j]
	}
barplot(mdata,names.arg = year, col =gray.colors(dat$na))
title(main='Summter Trawl Survey')

mdata <- matrix(0,dat$na,ny)
# Create Observer freq figure 
bdata <- bdata <- ono+ooo 
for (j in 1:dat$nyo) {
	mdata[,dat$io[j]-dat$fyear+1] <- bdata[,j]
	}
barplot(mdata,names.arg = year, col =gray.colors(dat$na))
title(main='Summer Observer Survey')

#mdata <- matrix(0,6,ny)
# Create Summer Pot length freq figure 
#bdata <- bdata <- dat$onp+dat$oop 
#for (j in 1:dat$nyp) {
#	mdata[,dat$ip[j]] <- bdata[,j]
#	}
#barplot(mdata,names.arg = year[-ny], col=color)
#legend('left',legend=c(6,5,4,3,2,1), fill=c('gray90','black','darkgray','dimgray','gray','white') , bg='white')
#title(main='Summter Pot Survey')
#mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)
plot.new()
legend('top',legend=len.c,fill=gray.colors(dat$na), bg='white', horiz=TRUE)


############################################################################
#  Figure for Length Frequency Plot
############################################################################
#par(mfrow=c(2,3),mar = c(2, 1.5, 1.5, 1.5),oma=c(4,1,2,1))
#lengthclass  <- c(74.0, 84.0, 94.0, 104, 114, 124,134)
# Create commercial length vs predicted length frequencies
#cdata <- rowSums((dat$onc+dat$ooc)*dat$tc)
#ccdata <- c(cumsum(cdata)/sum(cdata),1)
#pdata <- rowSums((rept$enc + rept$eoc)*dat$tc)
#cpdata <- c(cumsum(pdata)/sum(pdata),1)
#plot(lengthclass, ccdata, type='s', main='Commercial Harvest', xlab = 'Length', ylab='')
#lines(lengthclass, cpdata, type='s',lty=5, col='red')

# Create Winter length vs predicted length freq figure 
#cdata <- rowSums((dat$onw+dat$oow)*dat$sw)
#ccdata <- c(cumsum(cdata)/sum(cdata),1)
#pdata <- rowSums((rept$enw + rept$eow)*dat$sw)
#cpdata <- c(cumsum(pdata)/sum(pdata),1)
#plot(lengthclass, ccdata, type='s', main='Winter Pot Survey', xlab = 'Length', ylab='')
#lines(lengthclass, cpdata, type='s', lty=5, col='red')

# Create Trawl length vs predicted length freq figure 
#cdata <- rowSums((dat$ont+dat$oot)*dat$tt)
#ccdata <- c(cumsum(cdata)/sum(cdata),1)
#pdata <- rowSums((rept$ent)*rept$ett)
#cpdata <- c(cumsum(pdata)/sum(pdata),1)
#plot(lengthclass, ccdata, type='s', main='Trawl Survey', xlab = 'Length', ylab='')
#lines(lengthclass, cpdata, type='s', lty=5, col='red')

# Create Observer survey length vs predicted length freq figure 
#cdata <- rowSums((dat$ono+dat$ooo)*dat$so)
#ccdata <- c(cumsum(cdata)/sum(cdata),1)
#pdata <- rowSums((rept$eno + rept$eoo)*dat$so)
#cpdata <- c(cumsum(pdata)/sum(pdata),1)
#plot(lengthclass, ccdata, type='s', main='Summer Observer', xlab = 'Length', ylab='')
#lines(lengthclass, cpdata, type='s', lty=5, col='red')

# Create Pot survey length vs predicted length freq figure 
#cdata <- rowSums((dat$onp+dat$oop)*dat$tp)
#ccdata <- c(cumsum(cdata)/sum(cdata),1)
#pdata <- rowSums((rept$enp + rept$eop)*rept$etp)
#cpdata <- c(cumsum(pdata)/sum(pdata),1)
#plot(lengthclass, ccdata, type='s', main='Summer Pot Survey', xlab = 'Length', ylab='')
#lines(lengthclass, cpdata, type='s', lty=5, col='red')
#mtext("Length frequency: Observed (black) vs. Predicted (red)", side=3,outer=T)
#dev.off()

############################################################################
#  Create Bubble Plot
############################################################################
par(mfrow=c(3,2), mar=c(2,1.2,1.2,1.2),oma=c(3,1,2,1),mgp=c(0,.25,0))
bubbleplot<-function(){
test2$radius <- .3*sqrt(abs(test2$fbdata)/pi) 
symbols(test2$classes,test2$years,circles=test2$radius, bg=ifelse(test2$fbdata<0,'white','black'),fg=ifelse(test2$fbdata==0,'white','black'), tck = -0.02, inches = FALSE,ylim=rev(range(test2$years)),xlab = '', ylab='')
}
# egr 1
bdata <- ptag1-rept$egr1
fbdata <- c(bdata[1,],bdata[2,],bdata[3,],bdata[4,],bdata[5,])
fbdata <- ifelse(fbdata==0,NA,fbdata)
classes <- rep(1:6,5)
years <- sort(rep(1:5,6))
test2 <- data.frame(cbind(years,classes,fbdata))
bubbleplot()
title(main='1980-92: Recovery after 1 year')
# egr 12
bdata <- ptag12-rept$egr12
fbdata <- c(bdata[1,],bdata[2,],bdata[3,],bdata[4,],bdata[5,])
fbdata <- ifelse(fbdata==0,NA,fbdata)
classes <- rep(1:6,5)
years <- sort(rep(1:5,6))
test2 <- data.frame(cbind(years,classes,fbdata))
bubbleplot()
title(main='1993-2014: Recovery after 1 year')
# egr 1
bdata <- ptag2-rept$egr2
fbdata <- c(bdata[1,],bdata[2,],bdata[3,],bdata[4,],bdata[5,])
fbdata <- ifelse(fbdata==0,NA,fbdata)
classes <- rep(1:6,5)
years <- sort(rep(1:5,6))
test2 <- data.frame(cbind(years,classes,fbdata))
bubbleplot()
title(main='1980-92: Recovery after 2 year2')
# egr 1
bdata <- ptag22-rept$egr22
fbdata <- c(bdata[1,],bdata[2,],bdata[3,],bdata[4,],bdata[5,])
fbdata <- ifelse(fbdata==0,NA,fbdata)
classes <- rep(1:6,5)
years <- sort(rep(1:5,6))
test2 <- data.frame(cbind(years,classes,fbdata))
bubbleplot()
title(main='1993-2014: Recovery after 2 years')
# egr 1
bdata <- ptag3-rept$egr3
fbdata <- c(bdata[1,],bdata[2,],bdata[3,],bdata[4,],bdata[5,])
fbdata <- ifelse(fbdata==0,NA,fbdata)
classes <- rep(1:6,5)
years <- sort(rep(1:5,6))
test2 <- data.frame(cbind(years,classes,fbdata))
bubbleplot()
title(main='1980-92: Recovery after 3 year')
# egr 1
bdata <- ptag32-rept$egr32
fbdata <- c(bdata[1,],bdata[2,],bdata[3,],bdata[4,],bdata[5,])
fbdata <- ifelse(fbdata==0,NA,fbdata)
classes <- rep(1:6,5)
years <- sort(rep(1:5,6))
test2 <- data.frame(cbind(years,classes,fbdata))
bubbleplot()
title(main='1993-2014: Recovery after 3 years')
mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", cex = 0.9, side=1,line=1, outer=T)



dev.off()



















