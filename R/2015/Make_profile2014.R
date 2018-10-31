###############################################################################################
#   Crab Data profile likelihood construction for 2013 
#   This program read original data, modify save, run ADMB get likelihood value, and repeat this.
###############################################################################################
library(PBSadmb)
admodeldir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2014/Model/')
#setwd(admodeldir)

###############################################################################################
#   Crab Data construction 
###############################################################################################
# Crab data reading routine 
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/read_crab_data2013.R')
# Crab ouptutdata reading routine 
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/read_rep_data2.R')

###############################################################################################
#   1.1 Define input file name and output file name and locations. 
###############################################################################################
# infun: input file name 
datadir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2014/profile/')
datname <- paste('ns3n2014R_M.dat')
# R-data file
datadir2 <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/')
infn <- paste(datadir2,'ns3n2014R.txt',sep='')

###############################################################################################
#   1.1.1 ADMB program file locations. 
###############################################################################################
# set model prefix: ADMB mode file prefix name
prefix <- paste(admodeldir,'ns3n2014_1_4',sep='')
# set retro data prefix

# Read crab data into R List file 
dat <- datatoR(infn)

########  Add weights to Summer pot survey likelihood ############################################
dat$SDRec <- 0.5
dat$log_initp  <- 8
dat$initp_phase <- 1
dat$qtno_phase <- 1
dat$M_phase <- -1
dat$ms_phase <- -1
dat$lamc <- 1
dat$rm_phase <- -1
dat$rm2 <- 1

 
###############################################################################################
#   1.2  set output 
###############################################################################################
# nrep: the number or replicate step
# M.prof: vecotr that store total likelihood 
# P.prof: matrix that store parameter varlue of each iteration 
# tf1.prof: matrix that store likelihood of each value

nrep <- 20
maxgrd <- rep(0,nrep)
M.prof <- rep(0,nrep)
P.prof <- matrix(0,ncol=25,nrow=nrep)
tf1.prof <- matrix(0,ncol=10,nrow=nrep)
MMB <- matrix(0,ncol= (dat$ny +1), nrow=nrep)
# Determine minimum (s1) value and increment (incr) of parameter changing  
s1 <- 0
incr <- 0.5
outdatname <- paste('profile_M_')
varname <- 'rm2'
###############################################################################################
#   2.0  Run profile likelihood analyses
###############################################################################################
for (i in 1:nrep){
# Determine parameters changing
dat$ms <- 1
dat$rm2 <- s1 + incr*i
#dat$SDRec <- s1 + incr*i

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
par1 <- read.table(st, sep="",header=T)
pars <- data.frame(t(par1[,3]))
names(pars) <- t(par1[,2])
pars <- pars[,!(names(pars) %in% c('log_relrec','flnp','legaln','legalb','mmb0'))]
st2 <- paste(datadir,outdatname,(s1 + incr*i),'.std',sep="")
file.copy(st, st2)
# Run data read routin
#rept <- readReport(fn)
test <- paste(prefix,'.par',sep="") 
test1 <- as.double(scan(test,nlines=1,quiet=T,what="")) 
maxgrad[i] <- test1[16]
MMB[i,] <- as.double(scan(fn,skip=86,nlines=1,nmax=(dat$ny +1),quiet=T,what="")) 
M.prof[i] <- as.double(scan(fn,skip=96,nlines=1,quiet=T,what="")) 
tf1.prof[i,] <- as.double(scan(fn,skip=98,nlines=1,quiet=T,what="")) 
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
M <- seq(s1+incr, s1+nrep*incr, by=incr)
plot(M, maxgrad,type='l', main='maximum gradient')
plot(M, M.prof,type='l', main='Total negative log likelihood')
plot(M, tf1.prof[,1],type='l', main='Trawl survey')
#plot(M, tf1.prof[,2],type='l', main='Pot survey')
plot(M, tf1.prof[,3],type='l', main='Commercial cpue')
plot(M, -tf1.prof[,4],type='l',main='Trawl length comp')
#plot(M, -tf1.prof[,5],type='l',main='Summer Pot size')
plot(M, -tf1.prof[,6],type='l',main='Winter Pot length comp')
plot(M, -tf1.prof[,7],type='l',main='Summer Commercial length comp')
plot(M, tf1.prof[,8],type='l', main='recruits')	
plot(M, -tf1.prof[,9],type='l',main='Observer length comp')

# Plot likelihood profile standardized 
par(mfrow=c(1,1),mar = c(5, 5, 5, 5))
plot(M, M.prof-min(M.prof,na.rm=TRUE),type='l', lwd=2, main='Total negative log likelihood',xlab = varname, ylab='Likelihood',ylim=c(0,8))
lines(M, tf1.prof[,1]-min(tf1.prof[,1],na.rm=TRUE),type='l')
#lines(M, tf1.prof[,2]-min(tf1.prof[,2]),lty=2)
lines(M, tf1.prof[,3]-min(tf1.prof[,3],na.rm=TRUE),lty=2)
lines(M, -tf1.prof[,4]-min(-tf1.prof[,4],na.rm=TRUE),lty=3,lwd=2)
#lines(M, -tf1.prof[,5]-min(-tf1.prof[,5]),lty=2,col=2)
lines(M, -tf1.prof[,6]-min(-tf1.prof[,6],na.rm=TRUE),lty=6,col=2)
lines(M, -tf1.prof[,7]-min(-tf1.prof[,7],na.rm=TRUE),lty=2,lwd=2,col=4)
lines(M, tf1.prof[,8]-min(tf1.prof[,8],na.rm=TRUE),lty=4,col=4)
lines(M, -tf1.prof[,9]-min(-tf1.prof[,9],na.rm=TRUE),lty=5,col=4)
legend("topright", bty ='n', lwd=c(2,1,1,2,1,2,1,1), lty=c(1,1,2,3,6,2,4,5),col=c(1,1,1,1,2,4,4,4),
legend = c('Total','Trawl survey','Commercial cpue','Trawl length comp','Winter Pot length comp','Summer Commercial length comp','recruits','Observer length comp') )


# Plot parameter changes
par(mfrow=c(5,5),mar = c(2, 2, 2, 2))
    for (j in 1:length(pars)) {
    plot(M, P.prof[,j],type='l',main=names(pars[j]),xlab = varname)
    }

# Plot parameter changes
year <- seq(1976,1975+dat$ny+1)
par(mfrow=c(5,5),mar = c(2, 2, 2, 2))
    for (j in 1:nrep) {
    plot(year, MMB[j,]/1000,type='l',main=paste(varname,'=',M[j]),ylim = c(2,15))
    bmsy <- mean(MMB[j,which(year==1980):(dat$ny+1)])
    abline(h=bmsy/1000, lwd=1, lty=1)
    }

# Plot likelihood profile standardized 
par(mfrow=c(1,1),mar = c(5, 5, 5, 5))
    plot(year, MMB[1,]/1000,type='l',ylab='MMB x 1000',ylim = c(0,15))
   for (j in 3:nrep) {
    lines(year, MMB[j,]/1000,type='l')
    }



		