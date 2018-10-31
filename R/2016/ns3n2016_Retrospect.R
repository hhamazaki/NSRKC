library(mgcv)
library(PBSadmb)

###############################################################################################
#   Crab Data construction 
###############################################################################################
rm(list=ls(all=TRUE))
admodeldir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2016/Model/')
setwd(admodeldir)
firstyear <- 1976
lastyear <- 2015
year <- seq(firstyear,lastyear)
# Define number of years 
ny <- length(year)
###############################################################################################
#   1.0 Read source file 
###############################################################################################
# Read read.R for program use
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/ADMB_rw_data.R')

source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Prediction_model_dataset/Length_summary_freq.r')

# Make retro data
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/2016/make_retro_data2016.R')
###############################################################################################
#   1.1 Define input data file name and output file name and locations. 
###############################################################################################
# infun: input file name 
datadir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/')
datadir2 <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2016/Model_Eval/')
infn <- paste(datadir,'ns3n2016Rtagc3.txt',sep='')
# outfn: output file name  
outfn <- paste(datadir2,'ns3n2016R',sep='')

###############################################################################################
#   1.1.1 ADMB program file locations. 
###############################################################################################
# ADMB model file name
prefix <- paste(admodeldir,'ns3n2016_Feb1c3',sep='')
#pdffile <- paste(admodeldir,'ns3n2014retrospect.pdf',sep='')

# fnpar: ADMB parameter ouput file name 
fn <- paste(prefix,'.rep',sep='')
fnpar <- paste(prefix,'.std',sep='')


###############################################################################################
#   1.2 Red ADMB data into R  
###############################################################################################
# Read crab data into R List file 
dat<-datatoR(infn)

###############################################################################################
#   1.3 Update data  
###############################################################################################
dat$na <- dim(lengthclass)[1]
dat$ra <- length(lengthclass[(which(lengthclass$lenclass < 94)),])
dat$sla <- length(lengthclass[(which(lengthclass$lenclass < 104)),])
dat$qyear <- 1992
dat$ntag1 <- dim(t93)[1]
dat$ntag2 <- dim(t94)[1]
dat$slm <- min.s + (inc.s-1)/2
dat$slt <- inc.s
dat$lg <- p.legal
dat$wm <- lwp
dat$msn <- c(length(lengthclass[(which(lengthclass$lenclass >= 124)),]),length(lengthclass[(which(lengthclass$lenclass < 84)),]) )
dat$ont <- as.matrix(trwl.n[,-c(1:2)])
dat$oot <- as.matrix(trwl.o[,-c(1:2)])
dat$onw <- as.matrix(wntr.n[,-c(1:2)])
dat$oow <- as.matrix(wntr.o[,-c(1:2)])
dat$onc <- as.matrix(comh.n[,-c(1:2)])
dat$ooc <- as.matrix(comh.o[,-c(1:2)])
dat$ono <- as.matrix(obs.n[,-c(1:2)])
dat$ooo <- as.matrix(obs.o[,-c(1:2)])
dat$tag_recov1 <- as.matrix(t93[,1:5])
dat$tag_recov2 <- as.matrix(t94[,1:5])
 
########  Add weights to Summer pot survey likelihood ############################################
dat$SDRec <- 0.5
dat$SDW <- 0.3
dat$log_initp  <- 8
dat$initp_phase <- 1
dat$qtno_phase <- 1
dat$M_phase <- -1
dat$ms <- 3.6
dat$ms_phase <- 1
dat$lamc <- 1
dat$lamw <- 0
dat$lawp <- 1
dat$latag <- 0.5
dat$st_est <- 1
dat$nst <- 1
dat$nsc <- 1

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
fn <- paste(prefix,'.rep',sep="")
runAD(prefix=prefix, argvec = argvec,logfile=TRUE, add = FALSE, verbose=FALSE)
# set output 
MMB[1:(ny+1-i),i+1] <- as.double(scan(fn,skip=97,nlines=1,quiet=T,what="",nmax=ny+1)) 
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

