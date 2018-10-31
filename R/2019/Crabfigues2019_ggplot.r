############################################################################
#  Figure producing R function rutine
############################################################################
rm(list=ls(all=TRUE))
library(PBSadmb)
library(mgcv)
library(officer)
library(ggplot2)
library(cowplot)
library(reshape2)


############################################################################
#   1.0 Read source file 
############################################################################
# Read read.R for program use
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/ADMB_rw_data.R')
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Prediction_model_dataset/Crab_model_data_2019.R')

# Crab Julyschedule change routine
#source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/make_July_data.R')

############################################################################
#   1.1 Define input file name and output file name and locations. 
############################################################################
# infun: input file name 
datadir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2019/Model/')
admodeldir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2019/Model/')
setwd(admodeldir)

############################################################################
#   Activate either, based on model you use 
############################################################################
#outfn <- paste(datadir,'ns3n2016R_5mm.dat',sep='')
outfn <- paste0(datadir,'ns3n2019_test.dat')
pdffile <- paste0(admodeldir,'ns3n2019.pdf')

############################################################################
#   1.1.1 ADMB program file locations: Choose prefix based on your modeling 
############################################################################
# ADMB model file name
prefix <- paste(admodeldir,'ns3n2019_Feb1_tr_retain',sep='')

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
dat$lyear <- 2018
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
# The number of winter Commercial retain survey 
dat$nywc <- length(com.w.yrs)
# The number of observer discards survey 
dat$nyod <- length(ob.d.yrs)
# The number of observer total catch survey 
dat$nyot <- length(ob.t.yrs)
# The number of spring tagging survey 
dat$nysp <- length(sp.yrs)
# Row number of tag data 
dat$ntag1 <- dim(t1w)[1]
# Row number of tag data 
dat$ntag2 <- dim(t2w)[1]
# summer commercial fishery event indices: 1)year large trawl fishery ended, 
# 2) escape mechanism installed,
# 3) year summer commrecial fishery CW>5 inch crab was accepted by buyers
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
# Winter Pots s length
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
# Observer discard survey years conducted 
dat$iod <- ob.d.yrs
# Summer Commercial Observer discards newshell length
dat$nonod <- as.matrix(ob.d.size[ob.d.size$Shell==1,-c(1:2)])
# Summer Commercial Observer discards oldshell length
dat$noood <- as.matrix(ob.d.size[ob.d.size$Shell==2,-c(1:2)])
# Observer total catch survey years conducted 
dat$iot <- ob.t.yrs
# Summer Commercial Observer total catch newshell length
dat$nonot <- as.matrix(ob.t.size[ob.t.size$Shell==1,-c(1:2)])
# Summer Commercial Observer total catch oldshell length
dat$nooot <- as.matrix(ob.t.size[ob.t.size$Shell==2,-c(1:2)])
# Winter Commercial survey years conducted 
dat$icw <- com.w.yrs
# Winter Commercial catch newshell length
dat$noncw <- as.matrix(com.size.w[com.size.w$Shell==1,-c(1:2)])
# Summer Commercial Observer total catch oldshell length
dat$noocw <- as.matrix(com.size.w[com.size.w$Shell==2,-c(1:2)])
# Spring Tagging survey years conducted  
dat$isp <- sp.yrs
# Spring Tagging survey newshell length
dat$nonsp <- as.matrix(sp.size.w[sp.size.w$Shell==1,-c(1:2)])
# Spring Tagging survey oldshell length
dat$noosp <- as.matrix(sp.size.w[sp.size.w$Shell==2,-c(1:2)])
# Tag recovery data first period
dat$tag_recov1 <- as.matrix(t1w)
# Tag recovery data second period
dat$tag_recov2 <- as.matrix(t2w)

########  Add weights to Summer pot survey likelihood ######################
# SD for recruitment and Winter CPUE
dat$SD <- c(0.5,0.3)

dat$base<- c(1,-1,-1,-1,-1,1)  
# 1: Summer Commercial logistic function: 1: one parameter, 2: two parameter
# 2: Estimate Summer Commercial retention curve -1: No, 1: Yes (phase)
# 3: Estimate Winter Commercial retention curve -1: No, 1: Yes (phase)
# 4: Estimate Anual Molting -1: No, 1: Yes (phase)
# 5: Estimate M  -1: No, 1: yes (phase)  
# 6: M of largest 2 classes: 0: separate, 1: same   
dat$nsel <- c(1,1,1,0,0)  
# Number of selectivity functions: Periods
# 1: Trawl survey
# 2: Summer Commercial
# 3: Summer Commercial Retention 
# 4: Spare space
# 5: Spare space 

dat$likew <- c(0.5,0,0,0,0)  
# Likelihood control 
# 1: Tagging likelihood weight
# 2: Sprint Tagging
# 3: Spare space 
# 4: Spare space
# 5: Spare space
 
dat$pwh <- 0.16
# Proportion of Winter Commercial catch 

#---------------------------------------------------------------------------
#   1.3 Update data  
#---------------------------------------------------------------------------
makeADdata(dat,outfn)

#---------------------------------------------------------------------------
#   2.0 Run ADMB Model
#--------------------------------------------------------------------------- 
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


#---------------------------------------------------------------------------
#   3.0 Model Output and Diagnoses 
#--------------------------------------------------------------------------- 

############################################################################
#  Calculate Sample numbers 
############################################################################
# Calculate total  
  st <- colSums(dat$nont)+colSums(dat$noot) 
  sw <- colSums(dat$nonw)+colSums(dat$noow)
  sod <- colSums(dat$nonod)+colSums(dat$noood)
  sot <- colSums(dat$nonot)+colSums(dat$nooot)  
  sc <- colSums(dat$nonc)+colSums(dat$nooc)
  scw <- colSums(dat$noncw)+colSums(dat$noocw)  
  ssp <- colSums(dat$nonsp)+colSums(dat$noosp)
  ont <- dat$nont/st[col(dat$nont)]
  oot <- dat$noot/st[col(dat$noot)]
  onw <- dat$nonw/sw[col(dat$nonw)]
  oow <- dat$noow/sw[col(dat$noow)]
  onc <- dat$nonc/sc[col(dat$nonc)]
  ooc <- dat$nooc/sc[col(dat$nooc)]
  onod <- dat$nonod/sod[col(dat$nonod)]
  oood <- dat$noood/sod[col(dat$noood)]
  onot <- dat$nonot/sot[col(dat$nonot)]
  ooot <- dat$nooot/sot[col(dat$nooot)]
  oncw <- dat$noncw/scw[col(dat$noncw)]
  oocw <- dat$noocw/scw[col(dat$noocw)]  
  onsp <- dat$nonsp/ssp[col(dat$nonsp)]
  oosp <- dat$noosp/ssp[col(dat$noosp)]

############################################################################
#  Organize Tag recovery data 
############################################################################
#---------------------------------------------------------------------------
#   Calculate tag recovery effective sample size 
#--------------------------------------------------------------------------- 
   
tag.recap<-function(na,col.rec,dat1,dat2){
 temp.matrix1 <- matrix(0,na,na)
 temp.matrix2 <- matrix(0,na,na)
 n1 <- dim(dat1)[1]
 n2 <- dim(dat2)[1]
 for(i in 1:n1) 
  {
  temp.matrix1[dat1[i,1],dat1[i,2]] <- dat1[i,col.rec]
  }
 for(i in 1:n2) 
  {
  temp.matrix2[dat2[i,1],dat2[i,2]] <- dat2[i,col.rec]
  }  
 outmatrix <- temp.matrix1 + temp.matrix2
 outmatrix
}

# Calcurate the number of tag recoveries matrix 
tagrecap1 <- tag.recap(dat$na,3,dat$tag_recov1,dat$tag_recov2)
tagrecap2 <- tag.recap(dat$na,4,dat$tag_recov1,dat$tag_recov2)
tagrecap3 <- tag.recap(dat$na,5,dat$tag_recov1,dat$tag_recov2)

# calculate input sample size (sum of each row frequency)
nrecap <- c(rowSums(tagrecap1),rowSums(tagrecap2), 
            rowSums(tagrecap3))

# Change to proportions 
tagrecap1 <- tagrecap1/rowSums(tagrecap1)
tagrecap2 <- tagrecap2/rowSums(tagrecap2)
tagrecap3 <- tagrecap3/rowSums(tagrecap3)

# replace NaN to 0 
tagrecap1[is.nan(tagrecap1)] <- 0
tagrecap2[is.nan(tagrecap2)] <- 0
tagrecap3[is.nan(tagrecap3)] <- 0

# calculate effective sample size 
esn1 <- rowSums(rept$egr1*(1-rept$egr1))/rowSums((rept$egr1-tagrecap1)^2)
esn2 <- rowSums(rept$egr2*(1-rept$egr2))/rowSums((rept$egr2-tagrecap2)^2)
esn3 <- rowSums(rept$egr3*(1-rept$egr3))/rowSums((rept$egr3-tagrecap3)^2)
enrecap <- c(esn1, esn2, esn3)
# replace NaN to 0 
enrecap[is.nan(enrecap)] <- 0

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
 
#===========================================================================
#   4.1  Figures of Molting probability and Selectivity  
#===========================================================================
g_step <- function(g,main,xlab,ylab){
 g+geom_step()+coord_cartesian(ylim=c(0,1))+
 labs(title=main,x=xlab,y=ylab)
}
 
 
#---------------------------------------------------------------------------
#   4.1.1  Plot Molting probability  
#--------------------------------------------------------------------------- 
molt <- c(rept$molp[,1],rept$molp[dat$na,1])
trawl <- c(rept$selt[1,],rept$selt[1,dat$na])
wpot <- c(rept$selw,rept$selw[dat$na])
scl <- c(rept$selc,rept$selc[dat$na])
scr <- c(rept$selr,rept$selr[dat$na])
labels <- c(molt='Molting Probability',trawl='Trawl Selectivity',wpot='Winter Pote Selectivity',
            scl='Commercial Pot Seleictivity',scr='Commercial Retain Probability')
df <- data.frame(cbind(length.c,molt,trawl,wpot,scl,scr))
df.l <- melt(df, vraiable.name='key',value.names ='prob',id.vars='length.c')
ggplot(df.l,aes(length.c,value))+geom_step()+
facet_wrap(~variable,ncol=3,nrow=3, labeller=labeller(variable=labels))+
coord_capped_cart(bottom='both', left='both')



#---------------------------------------------------------------------------
#   4.1.2  Plot Trawl Selectivity   
#--------------------------------------------------------------------------- 
df <- data.frame(cbind(length.c,trawlpn))
g <- ggplot(df,aes(length.c, trawlpn))
g2 <- g_step(g,'Trawl Selectivity','','')

if(dat$nsel[1]==2) title(main='NOAA Trawl Selectivity') else title(main='Trawl Selectivity')
# ADFG or 2 trawl selectivity
trawlpa <- c(rept$selt[2,],rept$selt[2,dat$na])
if(dat$nsel[1]==2) plot(length.c, trawlpa,type='s', main='ADFG Trawl Selectivity', ylim=c(0,1))
#---------------------------------------------------------------------------
#   4.1.3  Winter Pot Selectivity   
#--------------------------------------------------------------------------- 
df <- data.frame(cbind(length.c,wotp))
g <- ggplot(df,aes(length.c, trawlpn))
g2 <- g_step(g,'Trawl Selectivity','','')




plot(length.c, wpotp,type='s', main='Winter pot Selectivity', ylim=c(0,1))
#---------------------------------------------------------------------------
#   4.1.4  Commercial Pot Selectivity   
#--------------------------------------------------------------------------- 
sc1p <- c(rept$selc,rept$selc[dat$na])
if(dat$nsel[2]==2)sc1p <- c(rept$selc[1,],rept$selc[1,dat$na])
plot(length.c, sc1p,type='s',ylim=c(0,1))
if(dat$nsel[2]==2) title(main='Commercial Selectivity 76-92') else title(main='Commercial Selectivity')
sc2p <- c(rept$selc1[2,],rept$selc1[2,dat$na])
if(dat$nsel[2]==2) plot(length.c, sc2p,type='s', main='Commercial Selectivity 93-', ylim=c(0,1))

#---------------------------------------------------------------------------
#   4.1.5  Commercial Retention Selectivity   
#--------------------------------------------------------------------------- 
sc1p <- c(rept$selr,rept$selr[dat$na])
if(dat$nsel[3]==2)sc1p <- c(rept$selr[1,],rept$selr[1,dat$na])
plot(length.c, sc1p,type='s',ylim=c(0,1))
if(dat$nsel[3]==2) title(main='Commercial Selectivity 76-92') else title(main='Commercial Retention ')
sc2p <- c(rept$selr[2,],rept$selr[2,dat$na])
if(dat$nsel[3]==2) plot(length.c, sc2p,type='s', main='Commercial Selectivity 93-', ylim=c(0,1))

plot_grid(,ncol=3,nrow=3)

#---------------------------------------------------------------------------
#   4.1.5  Plot Annual Molting probability  (Optional) 
#--------------------------------------------------------------------------- 
if(dat$rmol_phase==1)
	{
	par(mfrow=c(6,7),mar = c(2, 1.5, 1.5, 1.5),oma=c(4,1,2,1))
	for(i in 1:ny){
	moltp <- moltp <- c(rept$molp[,i],rept$molp[dat$na,i])
	plot(length.c, moltp,type='s', ylim=c(0,1),,main = year[i])
	}
	mtext("Annual Molting probability", side=3, outer=T)
	}

#===========================================================================
#   PLot QQ-plot
#===========================================================================
############################################################################
#  Calculate RMSE and Residual plots
############################################################################
trawl1 <- ifelse(is.null(pars$qtno),1,pars$qtno)*(rept$ett[1:6])
trawl2 <- ifelse(is.null(pars$qtad),1,pars$qtad)*(rept$ett[7:dat$nyt])
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
#  Plot Effective sample size 
############################################################################
par(mfrow=c(5,3),mar = c(3,3,2, 2),oma=c(3,1,2,1))

plot.efn<-function(estp,obp,sn,efn,maxs,tmain,year){
tefn <-colSums(estp*(1-estp))/colSums((estp-obp)^2)
tefnn <- ifelse(sn*efn>maxs,maxs,sn*efn)
# plot Implied Frequency
hist(tefn,breaks=10, main=tmain)
title(xlab = 'Implied', line = 2)
title(ylab = 'Frequency', line = 2)
abline(v=mean(tefn),lwd=2)
# plot Input vs. Implied
plot(tefnn,tefn, main=tmain, xlim = c(0,maxs),ylim=c(0,max(tefn)))
abline(lm(tefn~ tefnn -1),lty=2)
abline(0,1)
title(xlab = 'Input', line = 2)
title(ylab = 'Implied', line = 2)
# Plot Implied vs. year
plot(year, tefn,main=tmain)
title(xlab = 'Year', line = 2)
title(ylab = 'Implied', line = 2)
}


# Trawl Survey Figure
tmain <- 'Trawl survey'
et <- rept$ent+rept$eot
ot <- ont+oot
plot.efn(et,ot,st,dat$efn[1],dat$maxs[1],tmain,dat$it)

# Commercial Catch Figure
tmain <- 'Commercial Catch'
ec <- rept$enc+rept$eoc
oc <- onc+ooc
plot.efn(ec[,-c(1,16)],oc[,-c(1,16)],sc[-c(1,16)],dat$efn[2],dat$maxs[2],tmain,year[-c(1,16)])

# Winter Pot Survey Figure
tmain <- 'Winter pot survey'
ew <- rept$enw+rept$eow
ow <- onw+oow
plot.efn(ew,ow,sw,dat$efn[2],dat$maxs[2],tmain,dat$iw)

if(dat$base[2] < 0) { 
# Observer Survey Figure
tmain <- 'Observer survey Discards'
eo <- rept$enod+rept$eood
oo <- onod+oood
plot.efn(eo,oo,sod,dat$efn[2],dat$maxs[2],tmain,dat$iod)
} else {
# Observer Survey Figure
tmain <- 'Observer survey Total'
eo <- rept$enot+rept$eoot
oo <- onot+ooot
plot.efn(eo,oo,sot,dat$efn[2],dat$maxs[2],tmain,dat$iot)
}

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

if(dat$likew[2] == 0) { 
# Spring Pot Survey Figure
tmain <- 'Spring Pot survey'
# Spring Pot survey effective sample size
esp <- rept$ensp+rept$eosp
osp <- onsp+oosp

plot.efn(esp,osp,ssp,dat$efn[2],dat$maxs[2],tmain,dat$isp)
}


############################################################################
#  Figure for trends.
############################################################################
#---------------------------------------------------------------------------
#   Plot Trawl Survey
#--------------------------------------------------------------------------- 
par(mfrow=c(1,1), mar=c(5,5,4,4))
dat$stcpue[dat$stcpue==0]<-NA
ciwu1 <- exp(log(dat$stcpue)+2*sqrt(log((dat$secpue/dat$stcpue)^2+1)))
ciwl1 <- exp(log(dat$stcpue)-2*sqrt(log((dat$secpue/dat$stcpue)^2+1)))
ciwu <- exp(log(dat$stcpue)+2*sqrt(log((dat$secpue/dat$stcpue)^2+1)+pars$advar))
ciwl <- exp(log(dat$stcpue)-2*sqrt(log((dat$secpue/dat$stcpue)^2+1)+pars$advar))
plot(year, dat$stcpue, ylab = 'CPUE',  xlab = 'Year', ylim=c(0,7), bty='l')
arrows(year,y0=ciwu,y1=ciwl,code=0,col='red')
arrows(year,y0=ciwu1,y1=ciwl1,code=3,length=0.025, angle =90)
points(year, dat$stcpue, pch=21, bg='white')
rept$ecpue[rept$ecpue==0]<-NA
lines(year,rept$ecpue,lwd=2)
legend("topright",lty =c(1,0), lwd = c(2,0), pch=c(NA,19), legend = c('Predicted','Observed'))

# Figure 2. Total, Recruites, Legals estimated abundance 
#: Add projected number 
#rept$nps[1:5,ny+1] <- rept$npw[1:5,ny+1]*exp(-0.417*dat$M)
#rept$nps[6,ny+1] <- rept$npw[6,ny+1]*exp(-0.417*dat$M*dat$ms)
#rept$ops[1:5,ny+1] <- rept$opw[1:5,ny+1]*exp(-0.417*dat$M)
#rept$ops[6,ny+1] <- rept$opw[6,ny+1]*exp(-0.417*dat$M*dat$ms)

#---------------------------------------------------------------------------
#   Plot Total, Legal, Recruites
#--------------------------------------------------------------------------- 
total.abundance <- colSums(rept$npw+rept$opw)
plot(year.1, total.abundance/1000, ylab = 'Abundance (million crabs)', 
           typ = 'l', lwd=2, xlab = 'Year',ylim= c(0,10),bty='l')
legal.abundance <- colSums(rept$npw[(dat$sla+1):dat$na,]+rept$opw[(dat$sla+1):dat$na,])
lines(year.1, legal.abundance/1000, lty=2)
recruit.abundance <- colSums(rept$npw[1:dat$ra,]+rept$opw[1:dat$ra,])
lines(year.1, recruit.abundance/1000, lty=3)
legend("topright",lty = c(1,2,3), legend = c('total','legal','recruits'))
points(year.1[ny+1],total.abundance[(ny+1)]/1000,pch=19,col=1,cex=0.8)
points(year.1[ny+1],legal.abundance[(ny+1)]/1000,pch=19,col=1,cex=0.8)
points(year.1[ny+1],recruit.abundance[(ny+1)]/1000,pch=19,col=1,cex=0.8)

title(main='Modeled crab abundance Feb 01')

#---------------------------------------------------------------------------
#   Plot MMB
#--------------------------------------------------------------------------- 
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
tier <- ifelse(cofl>=1,'a','b')
# calculate Bmsy adjusted coefficient
fl <- ifelse(cofl>1,1,ifelse(0.25<cofl&cofl<=1,(cofl-0.1)/0.9,0))
# calculate FOF:
#rept$M <- 0.18

FOFL <- 0.18*fl
# Calculate Legal biomass = projected winter*fishery selectivity
#                           *prop legal*weight 
Legal.b <- with(rept,(npw[,p.y]+opw[,p.y])
			*dat$lg*rept$selcy[,p.y-1]*dat$wm)
Legal.s.b <- with(rept,npw[,p.y]+opw[,p.y])*(1-dat$lg)*rept$selcy[,p.y-1]*dat$wm			
# Calculate adjusted FOFL 
par1 <- 1-exp(-(FOFL +0.42*rept$M))
par2 <- 1-exp(-0.42*rept$M)
pw <- dat$pwh
FOFL.a <- par1 - par2*(1-pw*par1)/(1-pw*par2)
OFL.r <- sum(Legal.b*FOFL.a)
OFL.nr <- sum(Legal.s.b*(FOFL.a)*dat$hm[1])
ABC <- OFL.r*0.8

#OFLrb <- (1-exp(-FOFL))*pars$last_legal/1000
title(main='MMB Feb 01')
tex <- c(paste('BMSY ',round(bmsy/1000,3),' mil.lb'),
         paste('MMB ',round(rept$MMB[ny+1]/1000,3),' mil.lb'),
		 paste0('Tier 4',tier), 
		 paste('Legal B ',round(sum(Legal.b)/1000,3),' mil.lb'), 
		 paste('OFL ',round(OFL.r/1000,3),' mil.lb'),
		 paste('ABC ',round(ABC/1000,3),' mil.lb'))
#tex <- c(paste('BMSY ',round(bmsy/1000,3),' mil.lb'),paste('MMB ',round(rept$MMB[length(rept$MMB)]/1000,3),' mil.lb'),paste('Legal B ',round(rept$Legalb[length(rept$Legalb)]/1000,3),' mil.lb'))
legend("topright",tex,bty='n',xjust=0) 


#---------------------------------------------------------------------------
#   Plot Commercial CPUE 
#--------------------------------------------------------------------------- 
dat$stcpue[dat$stcpue==0]<-NA

ciwu1 <- exp(log(dat$stcpue)+2*sqrt(log((dat$secpue/dat$stcpue)^2+1)))
ciwl1 <- exp(log(dat$stcpue)-2*sqrt(log((dat$secpue/dat$stcpue)^2+1)))
ciwu <- exp(log(dat$stcpue)+2*sqrt(log((dat$secpue/dat$stcpue)^2+1)+pars$advar))
ciwl <- exp(log(dat$stcpue)-2*sqrt(log((dat$secpue/dat$stcpue)^2+1)+pars$advar))
ciw[ciw==0]<-NA
plot(year, dat$stcpue, ylab = 'CPUE',  xlab = 'Year', ylim=c(0,7), bty='l')
arrows(year,y0=ciwu,y1=ciwl,code=0,col='red')
arrows(year,y0=ciwu1,y1=ciwl1,code=3,length=0.025, angle =90)
points(year, dat$stcpue, pch=21, bg='white')
rept$ecpue[rept$ecpue==0]<-NA
lines(year,rept$ecpue,lwd=2)
legend("topright",lty =c(1,0), lwd = c(2,0), pch=c(NA,19), legend = c('Predicted','Observed'))
title(main='Summer commercial standardized cpue')
 
#---------------------------------------------------------------------------
#   Plot Catch data 
#--------------------------------------------------------------------------- 
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
#---------------------------------------------------------------------------
#   Commercial Length: predicted vs. Observed
#--------------------------------------------------------------------------- 
#par(mfrow=c(6,7),mar = c(0, 0, 0, 0),oma=c(4,4,4,4))
#for (i in c(2:15,17:length(year))){
#plot(len.c,onc[,i]+ooc[,i],ylim=c(0,0.7),pch=19, cex=0.8,tck = -0.05, 
#     axes=FALSE, bty='l')
#lines(len.c,rept$enc[,i]+rept$eoc[,i],lty=2)
#text(len.c[3],0.6,year[i],cex=1.0)
#if (i %in% c(2,9,17,24,31,38)) axis(2,at=seq(0,0.7,0.2))
#if (i %in% c(36:length(year))) axis(1,len.c)
#box()
#}
#mtext("commercial harvest length: observed vs predicted", side=3, outer=T)
#mtext("CL mm", side=1,line=2., outer=T)
#mtext("Proportion", side=2,line=2., outer=T)

#---------------------------------------------------------------------------
#   Commercial Length: predicted vs. Observed: New Shell & Old shell
#--------------------------------------------------------------------------- 
par(mfrow=c(6,7),mar = c(0, 0, 0, 0),oma=c(4,4,4,4))
for (i in c(2:15,17:length(year))){
plot(len.c,onc[,i],ylim=c(0,0.7),pch=19, cex=0.8,tck = -0.05, 
     axes=FALSE, bty='l')
lines(len.c,rept$enc[,i],lty=2)
points(len.c,ooc[,i],pch=19,col=2)
lines(len.c,rept$eoc[,i],lty=2,col=2)
text(len.c[3],0.6,year[i],cex=1.0)
if (i %in% c(2,9,17,24,31,38)) axis(2,at=seq(0,0.7,0.2))
if (i %in% c(37:length(year))) axis(1,len.c)
box()
}
mtext("commercial harvest length New Shell: observed vs predicted", side=3, outer=T)
mtext("CL mm", side=1,line=2., outer=T)
mtext("Proportion", side=2,line=2., outer=T)


#---------------------------------------------------------------------------
#   Winter & Spring Pot surve  Length: predicted vs. Observed
#--------------------------------------------------------------------------- 
par(mfrow=c(6,7),mar = c(0, 0, 0, 0),oma=c(4,4,4,4))
for (i in 1:dat$nyw){
plot(len.c,onw[,i]+oow[,i],ylim=c(0,0.5), pch=19, cex=0.8,tck = -0.05,
      axes=FALSE, bty='l', xlab = " ", ylab=" ")
lines(len.c,rept$enw[,i]+rept$eow[,i],lty=2)
text(len.c[3],0.45,dat$iw[i],cex=1.0)
if (i %in% c(1,8,15,22)) axis(2,at=seq(0,0.7,0.2))
if (i %in% c(21:dat$nyw)) axis(1,len.c)
box()
}
mtext("Winter pot length: observed vs predicted", side=3, outer=T)

#---------------------------------------------------------------------------
#   Winter & Spring Pot survey: New Shell  Length: predicted vs. Observed
#--------------------------------------------------------------------------- 
par(mfrow=c(6,7),mar = c(0, 0, 0, 0),oma=c(4,4,4,4))
for (i in 1:dat$nyw){
plot(len.c,onw[,i],ylim=c(0,0.5), pch=19, cex=0.8,tck = -0.05,
      axes=FALSE, bty='l', xlab = " ", ylab=" ")
lines(len.c,rept$enw[,i],lty=2)
points(len.c,oow[,i],pch=19, col=2)
lines(len.c,rept$eow[,i],lty=2,col=2)
text(len.c[3],0.45,dat$iw[i],cex=1.0)
if (i %in% c(1,8,15,22)) axis(2,at=seq(0,0.7,0.2))
if (i %in% c(21:dat$nyw)) axis(1,len.c)
box()
}
mtext("Winter pot length New Shell: observed vs predicted", side=3, outer=T)

#---------------------------------------------------------------------------
#   Winter & Spring Pot survey Old Shell Length: predicted vs. Observed
#--------------------------------------------------------------------------- 
par(mfrow=c(6,7),mar = c(0, 0, 0, 0),oma=c(4,4,4,4))
for (i in 1:dat$nyw){
plot(len.c,oow[,i],ylim=c(0,0.5), pch=19, cex=0.8,tck = -0.05,
      axes=FALSE, bty='l', xlab = " ", ylab=" ")
lines(len.c,rept$eow[,i],lty=2)
text(len.c[3],0.45,dat$iw[i],cex=1.0)
if (i %in% c(1,8,15,22)) axis(2,at=seq(0,0.7,0.2))
if (i %in% c(21:dat$nyw)) axis(1,len.c)
box()
}
mtext("Winter pot length Old Shell: observed vs predicted", side=3, outer=T)

#---------------------------------------------------------------------------
#   Trawl  Length: predicted vs. Observed
#--------------------------------------------------------------------------- 
par(mfrow=c(6,7),mar = c(0, 0, 0, 0),oma=c(4,4,4,4))
# Create Spring Pot survey length vs predicted length freq figure 
for (i in 1:dat$nysp){
x <- plot(len.c,onsp[,i]+oosp[,i],ylim=(c(0,0.5)), pch=19, cex=0.8,
     axes=FALSE,tck = -0.05, bty='l', xlab = " ", ylab=" ")
lines(len.c,rept$ensp[,i]+rept$eosp[,i],lty=2)
text(len.c[3],0.5,dat$isp[i],cex=1.0)
if (i %in% c(1)) axis(2,at=seq(0,0.7,0.2))
if (i %in% c(1:dat$nysp)) axis(1,len.c)
box()
}
mtext("Spring Pot survey: observed vs predicted", side=3, line = -46, outer=T)
mtext("CL mm", side=1,line=2., outer=T)
mtext("Proportion", side=2,line=2., outer=T)

#---------------------------------------------------------------------------
#   Trawl  Length: predicted vs. Observed
#--------------------------------------------------------------------------- 
par(mfrow=c(6,7),mar = c(0, 0, 0, 0),oma=c(4,4,4,4))
for (i in 1:dat$nyt){
x <- plot(len.c,ont[,i]+oot[,i],ylim=(c(0,0.5)), pch=19, cex=0.8,tck = -0.05,
        axes=FALSE, bty='l', xlab = " ", ylab=" ")
lines(len.c,rept$ent[,i]+rept$eot[,i],lty=2)
text(len.c[3],0.45,paste(dat$it[i],trawl$Agent[i]),cex=1.0)
#text(len.c[3],0.48,paste(trawl$Agent[i]),cex=1.0)
if (i %in% c(1,8,15)) axis(2,at=seq(0,0.5,0.2))
if (i %in% c(10:dat$nyt)) axis(1,len.c)
box()
}
mtext("Trawl length: observed vs predicted", side=3, outer=T)

#---------------------------------------------------------------------------
#   Trawl  Length: predicted vs. Observed New Shell
#--------------------------------------------------------------------------- 
par(mfrow=c(6,7),mar = c(0, 0, 0, 0),oma=c(4,4,4,4))
for (i in 1:dat$nyt){
plot(len.c,ont[,i],ylim=(c(0,0.5)), pch=19, cex=0.8,tck = -0.05,
        axes=FALSE, bty='l', xlab = " ", ylab=" ")
lines(len.c,rept$ent[,i],lty=2)
points(len.c,oot[,i],pch=19, col=2)
lines(len.c,rept$eot[,i],lty=2,col=2)

text(len.c[3],0.45,paste(dat$it[i],trawl$Agent[i]),cex=1.0)
#text(len.c[3],0.48,paste(trawl$Agent[i]),cex=1.0)
if (i %in% c(1,8,15)) axis(2,at=seq(0,0.5,0.2))
if (i %in% c(11:dat$nyt)) axis(1,len.c)
box()
}
mtext("Trawl length New Shell  & Old Shell: observed vs predicted", side=3, outer=T)
#---------------------------------------------------------------------------
#   Trawl  Length: predicted vs. Observed
#--------------------------------------------------------------------------- 
par(mfrow=c(6,7),mar = c(0, 0, 0, 0),oma=c(4,4,4,4))
for (i in 1:dat$nyt){
plot(len.c,oot[,i],ylim=(c(0,0.5)), pch=19, cex=0.8,tck = -0.05,
        axes=FALSE, bty='l', xlab = " ", ylab=" ")
lines(len.c,rept$eot[,i],lty=2)
text(len.c[3],0.45,paste(dat$it[i],trawl$Agent[i]),cex=1.0)
#text(len.c[3],0.48,paste(trawl$Agent[i]),cex=1.0)
if (i %in% c(1,8,15)) axis(2,at=seq(0,0.5,0.2))
if (i %in% c(10:dat$nyt)) axis(1,len.c)
box()
}
mtext("Trawl length Old Shell: observed vs predicted", side=3, outer=T)





#---------------------------------------------------------------------------
#   Observer  Length: predicted vs. Observed
#--------------------------------------------------------------------------- 
par(mfrow=c(6,7),mar = c(0, 0, 0, 0),oma=c(4,4,4,4))
for (i in 1:dat$nyo){
x <- plot(len.c,ono[,i]+ooo[,i],ylim=(c(0,0.6)), pch=19, cex=0.8,tck = -0.05, 
          axes=FALSE,bty='l', xlab = " ", ylab=" ")
lines(len.c,rept$eno[,i]+rept$eoo[,i],lty=2)
text(len.c[3],0.55,dat$io[i],cex=1.0)
if (i %in% c(1,8,15)) axis(2,at=seq(0,0.6,0.2))
if (i %in% c(6:dat$nyo)) axis(1,len.c)
box()
}
plot.new()
plot.new()
#---------------------------------------------------------------------------
#   Observer  Length: New Shell: predicted vs. Observed
#--------------------------------------------------------------------------- 
par(mfrow=c(6,7),mar = c(0, 0, 0, 0),oma=c(4,4,4,4))
for (i in 1:dat$nyod){
plot(len.c,onod[,i],ylim=(c(0,0.6)), pch=19, cex=0.8,tck = -0.05, 
          axes=FALSE,bty='l', xlab = " ", ylab=" ")
lines(len.c,rept$enod[,i],lty=2)
points(len.c,oood[,i], pch=19, col=2)
lines(len.c,rept$eood[,i],lty=2,col=2)

text(len.c[3],0.55,dat$iod[i],cex=1.0)
if (i %in% c(1,8,15)) axis(2,at=seq(0,0.6,0.2))
if (i %in% c(7:dat$nyod)) axis(1,len.c)
box()
}


#---------------------------------------------------------------------------
#   Observer  Length: New Shell: predicted vs. Observed
#--------------------------------------------------------------------------- 
par(mfrow=c(6,7),mar = c(0, 0, 0, 0),oma=c(4,4,4,4))
for (i in 1:dat$nyot){
plot(len.c,onot[,i],ylim=(c(0,0.6)), pch=19, cex=0.8,tck = -0.05, 
          axes=FALSE,bty='l', xlab = " ", ylab=" ")
lines(len.c,rept$enot[,i],lty=2)
points(len.c,ooot[,i], pch=19, col=2)
lines(len.c,rept$eoot[,i],lty=2,col=2)

text(len.c[3],0.55,dat$iot[i],cex=1.0)
#if (i %in% c(1,8,15)) axis(2,at=seq(0,0.6,0.2))
#if (i %in% c(7:dat$nyot) axis(1,len.c))
box()
}

#---------------------------------------------------------------------------
#   Observer  Length: Old Shell: predicted vs. Observed
#--------------------------------------------------------------------------- 
par(mfrow=c(6,7),mar = c(0, 0, 0, 0),oma=c(4,4,4,4))
for (i in 1:dat$nyo){
plot(len.c,ooo[,i],ylim=(c(0,0.6)), pch=19, cex=0.8,tck = -0.05, 
          axes=FALSE,bty='l', xlab = " ", ylab=" ")
lines(len.c,rept$eoo[,i],lty=2)
text(len.c[3],0.55,dat$io[i],cex=1.0)
if (i %in% c(1,8,15)) axis(2,at=seq(0,0.6,0.2))
if (i %in% c(6:dat$nyo)) axis(1,len.c)
box()
}

mtext("Discards length: observed vs predicted", side=3, line = -36, outer=T)
mtext("CL mm", side=1,line=2., outer=T)
mtext("Proportion", side=2,line=2., outer=T)

#---------------------------------------------------------------------------
#   Tag length  Length: predicted vs. Observed
#--------------------------------------------------------------------------- 
       
par(mfrow=c(6,(dat$na-1)/1),mar = c(0, 0, 0, 0),oma=c(4,4,4,4))
for (i in 1:(dat$na-1)){
x <- plot(len.c,tagrecap1[i,],ylim=c(0,1),pch=19, cex=0.8, tck = -0.05, 
		axes=FALSE,bty='l',ylab = '', xlab = '')
lines(len.c,rept$egr1[i,],lty=2)
legend('topleft',paste(len.c[i],'mm'),bty='n',cex=1.0)
if (i %in% c(1)) axis(2,at=seq(0,1,0.2))
axis(1,len.c)
box()
}
mtext("Tag recovery data observed vs predicted", cex = 0.9,side=3, line = 2, outer=T)
mtext("Recovery after 1 year", cex = 0.9, side=3, line = 0.5, outer=TRUE)
plot.new()
plot.new()
plot.new()
plot.new()
plot.new()
plot.new()
plot.new()
for (i in 1:(dat$na-1)){
x <- plot(len.c,tagrecap2[i,],ylim=c(0,1),pch=19, cex=0.8, tck = -0.05, 
	axes=FALSE,bty='l',ylab = '', xlab = '')
lines(len.c,rept$egr2[i,],lty=2)
legend('topleft',paste(len.c[i],'mm'),bty='n',cex=1.0)
if (i %in% c(1)) axis(2,at=seq(0,1,0.2))
axis(1,len.c)
box()
}
mtext("Recovery after 2 years", cex = 0.9, side=3, line = -18.5,outer=TRUE)
plot.new()
plot.new()
plot.new()
plot.new()
plot.new()
plot.new()
plot.new()
for (i in 1:(dat$na-1)){
x <- plot(len.c,tagrecap3[i,],ylim=c(0,1),pch=19, cex=0.8, tck = -0.05, 
	axes=FALSE,bty='l',ylab = '', xlab = '')
lines(len.c,rept$egr3[i,],lty=2)
legend('topleft',paste(len.c[i],'mm'),bty='n',cex=1.0)
if (i %in% c(1)) axis(2,at=seq(0,1,0.2))
axis(1,len.c)
box()
}
mtext("Recovery after 3 years", cex = 0.9, side=3, line = -37.5,outer=TRUE)


############################################################################
#  Figure for bubble plot
############################################################################
bubbleplot<-function(data1){
data1$radius <- 1.5*sqrt(abs(data1$fbdata)/pi) 
with(data1,symbols(years,classes,circles=radius, 
  bg=ifelse(fbdata<0,'white','black'), 
  inches = FALSE, xlim=c(dat$fyear,dat$lyear),xlab = " ", ylab=" "))
}

b.data <- function(omt1,omt2,emt1,emt2,yi,na,len.c){
bdata <- omt1+omt2 - (emt1 + emt2) 
fbdata <- as.vector(bdata)
ny <- length(yi)
classes <- rep(len.c,ny)
years <-1
for (i in 1:ny){
years <- c(years,rep(yi[i],na))
}
years <- years[-1]
data1 <- data.frame(cbind(years,classes,fbdata))
data1
}
par(mfrow=c(4,1), mar=c(2,2,2,1),oma=c(3,1,2,1))

#---------------------------------------------------------------------------
#   Commercial Harvest 
#--------------------------------------------------------------------------- 
data1 <- b.data(onc,ooc,rept$enc,rept$eoc,year,dat$na,len.c)
bubbleplot(data1[with(data1,years != 1976 &years != 1991),])
title(main='Commercial Harvest')

#---------------------------------------------------------------------------
#   Commercial Discards 
#--------------------------------------------------------------------------- 
data1 <- b.data(onod,oood,rept$enod,rept$eood,dat$iod,dat$na,len.c)
bubbleplot(data1)
title(main='Observer Survey')

#---------------------------------------------------------------------------
#   Trawl Survey NOAA & ADFG 
#--------------------------------------------------------------------------- 
data1 <- b.data(ont,oot,rept$ent,rept$eot,dat$it,dat$na,len.c)
# Add survey identifir
agent <-c(rep(1,6),rep(2,5),1,rep(2,3),1,2)
survey <-1
for (i in 1:dat$nyt){
survey <- c(survey,rep(agent[i],dat$na))
}
survey <- survey[-1]
t.all <- data.frame(cbind(data1,survey))
bubbleplot(t.all[t.all$survey==1,])
title(main='Trawl Survey NOAA')
bubbleplot(t.all[t.all$survey==2,])
title(main='Trawl Survey ADFG')

#---------------------------------------------------------------------------
#   Winter Pot Survey 
#--------------------------------------------------------------------------- 
data1 <- b.data(onw,oow,rept$enw,rept$eow,dat$iw,dat$na,len.c)
bubbleplot(data1)
title(main='Winter Pot Survey')

#---------------------------------------------------------------------------
#   Spring Pot Survey 
#--------------------------------------------------------------------------- 
data1 <- b.data(onsp,oosp,rept$ensp,rept$eosp,dat$isp,dat$na,len.c)
bubbleplot(data1)
title(main='Spring Pot Survey')



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



















