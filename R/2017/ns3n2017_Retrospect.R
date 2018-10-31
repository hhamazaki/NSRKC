library(mgcv)
library(PBSadmb)

###############################################################################################
#   Crab Data construction 
###############################################################################################
rm(list=ls(all=TRUE))
admodeldir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2017/Model/')
setwd(admodeldir)
firstyear <- 1976
lastyear <- 2016
year <- seq(firstyear,lastyear)
# Define number of years 
ny <- length(year)
###############################################################################################
#   1.0 Read source file 
###############################################################################################
# Read read.R for program use
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/ADMB_rw_data.R')
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Prediction_model_dataset/Crab_model_data.R')

# Make retro data
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/2017/make_retro_data2017.R')
###############################################################################################
#   1.1 Define input data file name and output file name and locations. 
###############################################################################################
# infun: input file name 
datadir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/')
datadir2 <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2017/Model/Model_Eval/')
infn <- paste(datadir,'ns3n2017new.dat',sep='')
# outfn: output file name  
outfn <- paste(datadir2,'ns3n2017R',sep='')

###############################################################################################
#   1.1.1 ADMB program file locations: Choose prefix based on your modeling 
###############################################################################################
# ADMB model file name
prefix <- paste(admodeldir,'ns3n2017_Feb1_tr',sep='')

# fnpar: ADMB parameter ouput file name 
fn <- paste(prefix,'.rep',sep='')
fnpar <- paste(prefix,'.std',sep='')


###############################################################################################
#   1.2 Red ADMB data into R  
###############################################################################################
# Read crab data into R List file 
dat<-list()
######## Create ADMB data ############################################
dat$fyear <- 1976
dat$lyear <- 2016
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
dat$ms6 <- 1.0
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

###############################################################################################
#   1.3 Update data  
###############################################################################################
########  Add weights to Summer pot survey likelihood ############################################
dat$SDRec <- 0.5
dat$SDW <- 0.3
dat$log_initp  <- 8
dat$initp_phase <- 1
dat$qtno_phase <- 1
dat$M_phase <- -1
dat$ms_phase <- 1
dat$rmol_phase <- -1
dat$lamc <- 1
dat$lamw <- 0
dat$lawp <- 1
dat$latag <- 0.5
dat$nst <- 1
dat$nsc <- 1
dat$smol <- -1
dat$ssc <- -1
dat$sst <- -1
dat$ssw <- -1
dat$sig <- -1
dat$ssth <- 3
dat$sthlike3 <- 0
dat$swm <- 1
dat$mol2p <- 1
dat$pwh <- 0.16

###############################################################################################
#   1.3 Create retrospective data  
###############################################################################################
# The file will be written from 0 to 
# Determine number of retrospective years
rn <- 10
# Run a model to prodcuce 0 to rn number of retrospective years
makeretrodata(dat,rn,outfn)

###############################################################################################
#   1.4  Run ADMB model and Extract output to R data 
###############################################################################################
# Make sure change data 
# Create empty data 
MMB <- matrix(0,ncol=rn+1, nrow=ny+1)
# Run the ADMB 
for (i in 0:rn){
argvec <- paste(paste(' -ind ',outfn,i,'.dat',sep=""),' -inp ',paste(prefix,'.par',sep=""),sep="")
#argvec <- paste(' -ind ',outfn,i,'.dat',sep="")
runAD(prefix=prefix, argvec = argvec,logfile=TRUE, add = FALSE, verbose=FALSE)
fn <- paste(prefix,'.rep',sep="")
# set output 
MMB[1:(ny+1-i),i+1] <- as.double(scan(fn,skip=121,nlines=1,quiet=T,what="",nmax=ny+1-i)) 
}
MMB

# Extract hind cast data for past years
hind <- MMB[(ny+1-rn):(ny),1]
hindpred <- numeric(rn)
for (i in rn:1){
hindpred[rn+1-i] <- MMB[ny+1-i,1+i]
}
hindpred


###############################################################################################
#   1.5  plot 
###############################################################################################
windows(record=TRUE)
#dev.copy(pdf, pdffile)
# Plot Rectrospecitve analysess
year <- seq(firstyear,lastyear+1)
plot(year, MMB[,1], col = 1, lwd = 3, type ='l', xaxt='n',ylab = 'MMB  ', main='Retrospective Analysis')
for (i in 1:rn){
a <- ny-i+1
lines(year[1:a],MMB[1:a,i+1], col = i+1 )
}
axis(side=1,seq(firstyear,lastyear+1,by=5)) 
hind.dif <- (hindpred-hind)/hind
Mohn.rho <- sum(hind.dif)
tex <- c(paste('Mohn rho ',round(Mohn.rho,3)))
legend("topright",tex,bty='n',xjust=0) 

#dev.off()
# Plot Retrospective vs. Predicted
plot(hind, hindpred, xlab='hindcast predicted legal (x1000)', ylab = 'retrospective MMB' )
abline(z <- lm(hind ~ hindpred -1))
coef(z)

minx <- min(min(hindpred),min(hind))
maxx <- max(max(hind)*1.2,max(hind))
year2 <- seq(lastyear-rn,lastyear-1)
plot(year2, hind, type = 'l', lty=1,ylim=c(minx,maxx), xlab='year', ylab = 'MMB (x1000)' )
lines(year2, hindpred, lty=2, col = 2)
legend('topright',c('hindcast MMB','predicted MMB'),lty=c(1,2),col=c(1,2),bty='n')
plot(hindpred, hind, xlab='predicted', ylab = 'hindcast MMB')

mean(hind.dif)
plot(year2, hind.dif)
write.csv(legal, paste(datadir,'legal.csv',sep=''))

