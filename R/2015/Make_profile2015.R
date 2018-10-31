###############################################################################################
#   Crab Data profile likelihood construction for 2013 
#   This program read original data, modify save, run ADMB get likelihood value, and repeat this.
###############################################################################################
library(PBSadmb)
admodeldir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2015/Model/')
#setwd(admodeldir)

###############################################################################################
#   Crab Data construction 
###############################################################################################
# Crab data reading routine 
# Read read.R for program use
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/read_crab_data.R')
# Crab ouptutdata reading routine 
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/read_rep_data2.R')
###############################################################################################
#   1.1 Define input file name and output file name and locations. 
###############################################################################################
# infun: input file name 
indatadir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/')
infn <- paste(indatadir,'ns3n2015Rtag.txt',sep='')
datadir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2015/profile/')
datname <- paste('ns3n2015R_M.dat')
###############################################################################################
#   1.1.1 ADMB program file locations. 
###############################################################################################
# ADMB model file name
prefix <- paste(admodeldir,'ns3n2015_Feb1_3',sep='')
# set retro data prefix

# Read crab data into R List file 
dat <- datatoR(infn)

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
dat$latag <- 0.5
#dat$efnw <- 0.1

 
###############################################################################################
#   1.2  set output 
###############################################################################################
# nrep: the number or replicate step
# M.prof: vecotr that store total likelihood 
# P.prof: matrix that store parameter varlue of each iteration 
# tf1.prof: matrix that store likelihood of each value

nrep <- 20
maxgrd <- rep(0,nrep)
x.axis <- rep(0,nrep)
M.prof <- rep(0,nrep)
P.prof <- matrix(0,ncol=27,nrow=nrep)
tf1.prof <- matrix(0,ncol=13,nrow=nrep)
MMB <- matrix(0,ncol= (dat$ny +1), nrow=nrep)
# Determine minimum (s1) value and increment (incr) of parameter changing  
s1 <- 1
incr <- 0.05
outdatname <- paste('profile_W_')
varname <- 'W'

###############################################################################################
#   2.0  Run profile likelihood analyses
###############################################################################################
for (i in 1:nrep){
# Determine parameters changing
#dat$ms <- 1.0
dat$latag <- s1 - incr*i
#dat$M <- s1 + incr*i
x.axis[i] <- s1 - incr*i
# Create a file with chaged parameter values
file.remove(paste(datadir,datname,sep=""))
out_file <- file(paste(datadir,datname,sep=""), open='a')
for (j in seq_along(dat)){
    write.table(paste("#",names(dat)[j]), file=out_file, sep="", dec=".", quote=FALSE, col.names=FALSE, row.names=FALSE)  
	write(t(dat[[j]]), file=out_file, sep=" ", 
	ncolumns= ifelse(is.numeric(dim(dat[[j]])[2]),dim(dat[[j]])[2],length(dat[[j]])))
	}
close(out_file) 

# Run the ADMB 
argvec <- paste(' -ind ',datadir, datname,' -inp ',paste(prefix,'.par',sep=""),sep="")
#argvec <- paste(' -ind ',datadir, datname,sep="")
runAD(prefix=prefix, argvec = argvec, verbose=FALSE)
fn <- paste(prefix,'.rep',sep="") 
fn2 <- paste(datadir,outdatname,(s1 + incr*i),'.rep',sep="")
file.copy(fn, fn2)
st <- paste(prefix,'.std',sep="") 
par1 <- read.table(st,header=T)
pars <- data.frame(t(par1[,3]))
names(pars) <- t(par1[,2])
pars <- pars[,!(names(pars) %in% c('log_relrec','flnp','legaln','legalb','mmb0'))]
st2 <- paste(datadir,outdatname,(s1 + incr*i),'.std',sep="")
file.copy(st, st2)
# Run data read routin
#rept <- readReport(fn)
#test <- paste(prefix,'.par',sep="") 
#test1 <- as.double(scan(test,nlines=1,quiet=T,what="")) 
#maxgrad[i] <- test1[16]
MMB[i,] <- as.double(scan(fn,skip=97,nlines=1,nmax=(dat$ny +1),quiet=T,what="")) 
M.prof[i] <- as.double(scan(fn,skip=107,nlines=1,quiet=T,what="")) 
tf1.prof[i,] <- as.double(scan(fn,skip=109,nlines=1,quiet=T,what="")) 
# convert data.frame to vector 
p <- as.vector(pars, mode='numeric')
# insert data to P.prof matrix
    for (j in 1:length(pars)){
    P.prof[i,j] <- p[j]
    }
}

###############################################################################################
#   3.0  Run profile likelihood analyses
###############################################################################################
windows(record = TRUE)

# Plot likelihood profile for total and each component
par(mfrow=c(4,3),mar = c(2, 2, 2, 2))
#plot(M, maxgrad,type='l', main='maximum gradient')
plot(x.axis, M.prof,type='l', main='Total negative log likelihood')
plot(x.axis, tf1.prof[,1],type='l', main='Trawl survey')
#plot(M, tf1.prof[,2],type='l', main='Pot survey')
plot(x.axis, tf1.prof[,3],type='l', main='Commercial cpue')
plot(x.axis, tf1.prof[,4]+tf1.prof[,5],type='l',main='Trawl length comp')
plot(x.axis, tf1.prof[,6]+tf1.prof[,7],type='l',main='Winter Pot length')
plot(x.axis, tf1.prof[,8]+tf1.prof[,9],type='l',main='Summer Commercial length comp')
plot(x.axis, tf1.prof[,10]+tf1.prof[,11],type='l', main='Observer length comp')	
plot(x.axis, tf1.prof[,12],type='l',main='Recruits')
plot(x.axis, tf1.prof[,13],type='l',main='Tag recovery')

# Plot likelihood profile standardized 
par(mfrow=c(1,1),mar = c(5, 5, 5, 5))
plot(x.axis, M.prof-min(M.prof,na.rm=TRUE),type='l', lwd=2, main='Total negative log likelihood', ylab='Likelihood',ylim=c(0,8))
lines(x.axis, tf1.prof[,1]-min(tf1.prof[,1],na.rm=TRUE),type='l')
#lines(M, tf1.prof[,2]-min(tf1.prof[,2]),lty=2)
lines(x.axis, tf1.prof[,3]-min(tf1.prof[,3],na.rm=TRUE),lty=2)
lines(x.axis, tf1.prof[,4]-min(-tf1.prof[,5],na.rm=TRUE),lty=3,lwd=2)
#lines(M, -tf1.prof[,5]-min(-tf1.prof[,5]),lty=2,col=2)
lines(x.axis, tf1.prof[,6]+tf1.prof[,7]-min(tf1.prof[,6]+tf1.prof[,7],na.rm=TRUE),lty=6,col=2)
lines(x.axis, tf1.prof[,9]+tf1.prof[,8]-min(tf1.prof[,9]+tf1.prof[,8],na.rm=TRUE),lty=2,lwd=2,col=4)
lines(x.axis, tf1.prof[,11]+tf1.prof[,10]-min(tf1.prof[,11]+tf1.prof[,10],na.rm=TRUE),lty=4,col=4)
lines(x.axis, tf1.prof[,12]-min(tf1.prof[,12],na.rm=TRUE),lty=5,col=4)
lines(x.axis, tf1.prof[,13]-min(tf1.prof[,13],na.rm=TRUE),lty=6,col=4)
legend("topright", bty ='n', lwd=c(2,1,1,2,1,2,1,1,1), lty=c(1,1,2,3,6,2,4,5,6),col=c(1,1,1,1,2,4,4,4,4),
legend = c('Total','Trawl survey','Commercial cpue','Trawl length comp','Winter Pot length comp','Summer Commercial length comp','recruits','Observer length comp','tag recovry') )


# Plot parameter changes
par(mfrow=c(5,5),mar = c(2, 2, 2, 2))
colnames(P.prof) <- names(pars)
    plot(x.axis, exp(P.prof[,'log_initpop']),type='l',main='Initial Population size')
    plot(x.axis, (P.prof[,'qtno']),type='l',main='Survey Q NOAA')
	
	

    for (j in 1:length(pars)) {
    plot(x.axis, (P.prof[,j]),type='l',main=names(pars[j]),xlab = varname)
    }
colnames(P.prof) <- names(pars)
	
# Plot likelihood profile standardized 
par(mfrow=c(1,1),mar = c(5, 5, 5, 5))
    plot(year, MMB[1,]/1000,type='l',ylab='MMB x 1000',ylim = c(0,15))
   for (j in 3:nrep) {
    lines(year, MMB[j,]/1000,type='l')
    }
	
##############################################################################################################
#  Plot Selectivity profile 
##############################################################################################################
par(mfrow=c(2,3),mar = c(2, 1.5, 1.5, 1.5),oma=c(4,1,2,1))	
sizeclass <- c(78.5, 88.5, 98.5, 108.5, 118.5, 128.5, 128.5)
lengthclass  <- c(74.0, 84.0, 94.0, 104, 114, 124, 134)
selc <-function(log.par){
sel <- 1/(1+exp(exp(log.par)*(128.5-sizeclass)+log(1.0/0.999-1.0)))
sel
}
p.selc <-function(log.par){
for (i in 2:nrep){
sel <- selc(log.par[i])
lines(lengthclass, sel,type='s') 
 }
} 
moltp <- 1-1/(1+exp(exp(P.prof[1,'log_mo1'])*(78.5-sizeclass)+log(1.0/0.001-1.0)))
plot(lengthclass, moltp,type='s', main='Molting Probability', ylim=c(0,1))
for (i in 2:nrep){
moltp <- 1-1/(1+exp(exp(P.prof[i,'log_mo1'])*(78.5-sizeclass)+log(1.0/0.001-1.0)))
lines(lengthclass, moltp,type='s') 
}
# Plot NOAA Trawl Survey 
trawlpn <- selc(P.prof[1,'log_st1'])
plot(lengthclass, trawlpn,type='s', main='NOAA Trawl Selectivity', ylim=c(0,1))
p.selc(P.prof[,'log_st1'])
# Plot ADFG Trawl Survey 
trawlpa <- selc(P.prof[1,'log_st3'])
plot(lengthclass, trawlpn,type='s', main='ADFG Trawl Selectivity', ylim=c(0,1))
p.selc(P.prof[,'log_st3'])
# Plot Winter Pot Survey 
wpotp <- selc(P.prof[1,'log_sw1'])
wpotp[6:7] <- (P.prof[1,'sw3'])
plot(lengthclass, wpotp,type='s', main='Winter pot Selectivity', ylim=c(0,1))
for (i in 2:nrep){
wpotp <- selc(P.prof[i,'log_sw1'])
wpotp[6:7] <- (P.prof[i,'sw3'])
lines(lengthclass, wpotp,type='s') 
}
# Plot Summer Commercial Catch 1976-1992 
sc1p <- selc(P.prof[1,'log_sc1'])
plot(lengthclass, sc1p,type='s', main='Commercial 77-92 Selectivity', ylim=c(0,1))
p.selc(P.prof[,'log_sc1'])

# Plot Summer Commercial Catch 1992-1 
sc3p <- selc(P.prof[1,'log_sc5'])
plot(lengthclass, sc1p,type='s', main='Commercial 93-13 Selectivity', ylim=c(0,1))
p.selc(P.prof[,'log_sc5'])

##############################################################################################################
#  Plot Selectivity profile by year 
##############################################################################################################
### Molting probability 
par(mfrow=c(5,5),mar = c(1.75,1.5,1.5,1.75),oma = c(3,3,3,3),cex=0.6) 
    for (i in 1:nrep) {
	moltp <- 1-1/(1+exp(exp(P.prof[1,'log_mo1'])*(78.5-sizeclass)+log(1.0/0.001-1.0)))
	plot(lengthclass, moltp,type='s', main=paste(varname,'=',x.axis[i]), ylim=c(0,1))
     }
mtext('Molting Probability', side = 3, line = 0, outer = TRUE)
mtext('Molting probability', side = 2, line = 1, outer = TRUE)
mtext("CL Length", side = 1, line = 1, outer = TRUE)

### NOAA Trawl Selectivity
par(mfrow=c(5,5),mar = c(1.75,1.5,1.5,1.75),oma = c(3,3,3,3),cex=0.6) 
    for (i in 1:nrep) {
    plot(lengthclass, selc(P.prof[i,'log_st1']),type='s',main=paste(varname,'=',x.axis[i]),ylim = c(0,1))
    }
mtext('NOAA Trawl Selectivity', side = 3, line = 0, outer = TRUE)
mtext('Selectivity', side = 2, line = 1, outer = TRUE)
mtext("CL Length", side = 1, line = 1, outer = TRUE)

### ADFG Trawl Selectivity
par(mfrow=c(5,5),mar = c(1.75,1.5,1.5,1.75),oma = c(3,3,3,3),cex=0.6) 
    for (i in 1:nrep) {
    plot(lengthclass, selc(P.prof[i,'log_st3']),type='s',main=paste(varname,'=',x.axis[i]),ylim = c(0,1))
    }
mtext('ADFG Trawl Selectivity', side = 3, line = 0, outer = TRUE)
mtext('Selectivity', side = 2, line = 1, outer = TRUE)
mtext("CL Length", side = 1, line = 1, outer = TRUE)

### Winter Pot Selectivity
par(mfrow=c(5,5),mar = c(1.75,1.5,1.5,1.75),oma = c(3,3,3,3),cex=0.6) 
    for (i in 1:nrep) {
	wpotp <- selc(P.prof[i,'log_sw1'])
	wpotp[6:7] <- (P.prof[i,'sw3'])
	plot(lengthclass, wpotp,type='s', main=paste(varname,'=',x.axis[i]), ylim=c(0,1))
     }
mtext('Winter Pot Selectivity', side = 3, line = 0, outer = TRUE)
mtext('Selectivity', side = 2, line = 1, outer = TRUE)
mtext("CL Length", side = 1, line = 1, outer = TRUE)

### Summer Commercial Catch Selectivity 1  
par(mfrow=c(5,5),mar = c(1.75,1.5,1.5,1.75),oma = c(3,3,3,3),cex=0.6) 
    for (i in 1:nrep) {
    plot(lengthclass, selc(P.prof[i,'log_sc1']),type='s',main=paste(varname,'=',x.axis[i]),ylim = c(0,1))
    }
mtext("Commercial 76-92 Selectivity", side = 3, line = 0, outer = TRUE)
mtext('Selectivity', side = 2, line = 1, outer = TRUE)
mtext("CL Length", side = 1, line = 1, outer = TRUE)
	
### Summer Commercial Catch Selectivity 2  
par(mfrow=c(5,5),mar = c(1.75,1.5,1.5,1.75),oma = c(3,3,3,3),cex=0.6) 
    for (i in 1:nrep) {
    plot(lengthclass, selc(P.prof[i,'log_sc5']),type='s',main=paste(varname,'=',x.axis[i]),ylim = c(0,1))
    }
mtext("Commercial 93-13 Selectivity", side = 3, line = 0, outer = TRUE)
mtext('Selectivity', side = 2, line = 1, outer = TRUE)
mtext("CL Length", side = 1, line = 1, outer = TRUE)




	
	
# Plot MMB changes
year <- seq(1976,1975+dat$ny+1)
par(mfrow=c(5,5),mar = c(2, 2, 2, 2))
    for (j in 1:nrep) {
    plot(year, MMB[j,]/1000,type='l',main=paste(varname,'=',x.axis[j]),ylim = c(2,15))
    bmsy <- mean(MMB[j,which(year==1980):(dat$ny+1)])
    abline(h=bmsy/1000, lwd=1, lty=1)
    }

# Plot likelihood profile standardized 
par(mfrow=c(1,1),mar = c(5, 5, 5, 5))
    plot(year, MMB[1,]/1000,type='l',ylab='MMB x 1000',ylim = c(0,15))
   for (j in 3:nrep) {
    lines(year, MMB[j,]/1000,type='l')
    }



