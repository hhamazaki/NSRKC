library(PBSadmb)
###############################################################################################
#   Crab Data construction 
###############################################################################################
rm(list=ls(all=TRUE))
admodeldir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2014/Model/')
setwd(admodeldir)
firstyear <- 1976
lastyear <- 2014
###############################################################################################
#   1.0 Read source file 
###############################################################################################
# Read read.R for program use
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/read_crab_data2013.R')
# Make retro data
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/2014/make_retro_data2014.R')

###############################################################################################
#   1.1 Define input data file name and output file name and locations. 
###############################################################################################
# infun: input file name 
datadir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/')
datadir2 <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2014/Model_Eval/')
infn <- paste(datadir,'ns3n2014Rtag2.txt',sep='')
# outfn: output file name  
outfn <- paste(datadir2,'ns3n2014R',sep='')


###############################################################################################
#   1.1.1 ADMB program file locations. 
###############################################################################################
# set moedel prefix
prefix <- paste(admodeldir,'ns3n2014_Jul1_c2',sep='')
pdffile <- paste(admodeldir,'ns3n2014retrospect.pdf',sep='')

###############################################################################################
#   1.2 Red ADMB data into R  
###############################################################################################
# Read crab data into R List file 
dat<-datatoR(infn)

###############################################################################################
#   1.2 Change data configuring into July schedule  
###############################################################################################
# Read crab data into R List file 
dat<-datatoR(infn)
# Winter commercial starts 1976 winter
dat$twc <- dat$twc[-1]
dat$tws <- dat$tws[-1]
dat$twst <- dat$twst[-1]
# Add zero to the last year of fishery
dat$twc <- c(dat$twc, dat$twc[dat$ny-1]) 
dat$tws <- c(dat$tws, dat$tws[dat$ny-1])
dat$twst <- c(dat$twst,dat$twst[dat$ny-1])

# Winter Pot survey is also shifted 1 year
dat$iw <- dat$iw-1

#dat$cv3[1:2] <- 0.01
#dat$ms <- 3.6
#dat$M <- 0.18
#dat$slt6 <- 1
#dat$ny <- 37

########  Add weights to Summer pot survey likelihood ############################################
dat$maxss <- 20
dat$maxsc <- 10
dat$SDRec <- 0.5
dat$SDW <- 0.3
dat$log_initp  <- 8
dat$initp_phase <- 1
dat$qtno_phase <- 1
dat$M_phase <- -1
dat$ms_phase <- -1
dat$lamc <- 1
dat$lamw <- 0
dat$lawp <- 1


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
MMB <- matrix(0,ncol=rn+1, nrow=dat$ny+1)
# Run the ADMB 
for (i in 0:rn){
#argvec <- paste(paste(' -ind ',outfn,i,'.dat',sep=""),' -inp ',paste(prefix,'.par',sep=""),sep="")
argvec <- paste(' -ind ',outfn,i,'.dat',sep="")
fn <- paste(prefix,'.rep',sep="")
runAD(prefix=prefix, argvec = argvec,logfile=TRUE, add = FALSE, verbose=FALSE)
# set output 
MMB[,i+1] <- as.double(scan(fn,skip=90,nlines=1,quiet=T,what="",nmax=dat$ny+1)) 
}
MMB

# Extract hind cast data for past years
hind <- MMB[(dat$ny+1-rn):(dat$ny),1]
hindpred <- numeric(rn)
for (i in rn:1){
hindpred[rn+1-i] <- MMB[dat$ny+1-i,1+i]
}
hindpred



###############################################################################################
#   1.5  plot 
###############################################################################################
windows(record=TRUE)
#dev.copy(pdf, pdffile)
# Plot Rectrospecitve analysess
year <- seq(firstyear,lastyear)
plot(year, MMB[,1], type ='l', xaxt='n',ylab = 'MMB  ', main='Retrospective Analysis')
for (i in 1:rn){
a <- dat$ny-i
lines(year[1:a],MMB[1:a,i+1], col = i)
}
axis(side=1,seq(firstyear,lastyear,by=5)) 
lines(seq(lastyear-rn,lastyear-1),hindpred, col = 'red', lwd = 2)
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
hind.dif <- (hindpred-hind)/hindpred
mean(hind.dif)
plot(year2, hind.dif)
write.csv(legal, paste(datadir,'legal.csv',sep=''))




