###############################################################################################
#   Crab Data profile likelihood construction for 2013 
#   This program read original data, modify save, run ADMB get likelihood value, and repeat this.
###############################################################################################
library(PBSadmb)
admodeldir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2013/Model/')
setwd(admodeldir)

###############################################################################################
#   Crab Data construction 
###############################################################################################
# Crab data reading routine 
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/read_crab_data2013.R')
# Crab ouptutdata reading routine 
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/read_rep_data.R')

###############################################################################################
#   1.1 Define input file name and output file name and locations. 
###############################################################################################
# infun: input file name 
datadir <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2014/profile/')
datname <- paste('ns3n2013R_M.dat')
# R-data file
datadir2 <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/')
infn <- paste(datadir2,'ns3n2012R.txt',sep='')

###############################################################################################
#   1.1.1 ADMB program file locations. 
###############################################################################################
# set model prefix: ADMB mode file prefix name
prefix <- paste(admodeldir,'ns3n2013cpue_rev2',sep='')
# set retro data prefix

# Read crab data into R List file 
dat <- datatoR(infn)

########  Add weights to Summer pot survey likelihood ############################################
dat$lamc <- 1
dat$log_initp  <- 8.0
dat$initp_phase <- 1
dat$maxss <- 20
dat$wsp <- c(0,0)
dat$qtno_phase <- 1
dat$qtad_phase <- -1 


###############################################################################################
#   set output 
###############################################################################################
nrep <- 12
M.prof <- rep(0,nrep)
P.prof <- matrix(0,ncol=25,nrow=nrep)
tf1.prof <- matrix(0,ncol=10,nrow=nrep)
s1 <- 0.18
incl <- 0.02
outdatname <- paste('profile_M_')
for (i in 1:nrep){
dat$ms <- 1
#dat$M <- 0.1
#dat$slt6 <- 1
dat$M <- s1 + incl*i

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
fn2 <- paste(datadir,outdatname,(s1 + incl*i),'.rep',sep="")
file.copy(fn, fn2)
st <- paste(prefix,'.std',sep="") 
st2 <- paste(datadir,outdatname,(s1 + incl*i),'.std',sep="")
file.copy(st, st2)
par1 <- read.table(st, sep="",header=T)
pars <- data.frame(t(par1[,3]))
names(pars) <- t(par1[,2])
# Run data read routin
#rept <- readReport(fn)
M.prof[i] <- as.double(scan(fn,skip=106,nlines=1,quiet=T,what="")) 
tf1.prof[i,] <- as.double(scan(fn,skip=108,nlines=1,quiet=T,what="")) 
P.prof[i,1] <- pars$log_q1
P.prof[i,2] <- pars$log_q2
#P.prof[i,3] <- pars$log_q3
P.prof[i,4] <- pars$log_initpop
P.prof[i,5] <- pars$log_recscale
P.prof[i,6] <- pars$r1
P.prof[i,7] <- pars$log_mo1
P.prof[i,8] <- pars$log_mo2
P.prof[i,9] <- pars$log_st1
P.prof[i,10] <- pars$log_st2
P.prof[i,11] <- pars$log_sw1
P.prof[i,12] <- pars$log_sw2
P.prof[i,13] <- pars$sw3
P.prof[i,14] <- pars$log_sc1
P.prof[i,15] <- pars$log_sc2
P.prof[i,16] <- pars$log_sc3
P.prof[i,17] <- pars$log_sc4
#P.prof[i,18] <- pars$log_sc5
#P.prof[i,19] <- pars$log_sc6
P.prof[i,20] <- pars$qtno
P.prof[i,21] <- pars$last_mmb
}
M.prof
#write.csv(M.prof, paste(datadir,'Mprog.csv',sep=''))
#tf1.prof
#write.csv(tf1.prof, paste(datadir,'tfprof.csv',sep=''))


par(mfrow=c(4,3),mar = c(2, 2, 2, 2))
M <- seq(s1+incl, s1+nrep*incl, by=incl)
plot(M, M.prof,type='l', main='Total negative log likelihood')
plot(M, tf1.prof[,1],type='l', main='Trawl survey')
#plot(M, tf1.prof[,2],type='l', main='Pot survey')
plot(M, tf1.prof[,3],type='l', main='Commercial cpue')
plot(M, -tf1.prof[,4],type='l',main='Trawl size')
#plot(M, -tf1.prof[,5],type='l',main='Summer Pot size')
plot(M, -tf1.prof[,6],type='l',main='Winter Pot size')
plot(M, -tf1.prof[,7],type='l',main='Summer Commercial size')
plot(M, tf1.prof[,8],type='l', main='recruits')	
plot(M, -tf1.prof[,9],type='l',main='Observer size')

par(mfrow=c(5,5),mar = c(2, 2, 2, 2))
plot(M, P.prof[,1],type='l',main='log_q1')
plot(M, P.prof[,2],type='l',main='log_q2')
plot(M, P.prof[,4],type='l',main='log_initpop')
plot(M, P.prof[,5],type='l',main='log_recscale')
plot(M, P.prof[,6],type='l',main='r1')
plot(M, P.prof[,7],type='l',main='log_mo1')
plot(M, P.prof[,8],type='l',main='log_mo2')
plot(M, P.prof[,9],type='l',main='log_st1')
plot(M, P.prof[,10],type='l',main='log_st2')
plot(M, P.prof[,11],type='l',main='log_sw1')
plot(M, P.prof[,12],type='l',main='log_sw2')
plot(M, P.prof[,13],type='l',main='sw3')
plot(M, P.prof[,14],type='l',main='log_sc1')
plot(M, P.prof[,15],type='l',main='log_sc2')
plot(M, P.prof[,16],type='l',main='log_sc3')
plot(M, P.prof[,17],type='l',main='log_sc4')
plot(M, P.prof[,20],type='l',main='qNOAA')
plot(M, P.prof[,21],type='l',main='last_MMB')



par(mfrow=c(1,1))
#M <- seq(s1+incl, s1+nrep*incl, by=incl)
#plot(M, M.prof-min(M.prof),type='l', main='Total negative log likelihood')
plot(M, tf1.prof[,1]-min(tf1.prof[,1]),type='l', col = 1, ylab='Likelihood',main='Total negative log likelihood')
#lines(M, tf1.prof[,2]-min(tf1.prof[,2]),lty=2)
lines(M, tf1.prof[,3]-min(tf1.prof[,3]),lty=5)
lines(M, -tf1.prof[,4]-min(-tf1.prof[,4]),lty=1,col=2)
#lines(M, -tf1.prof[,5]-min(-tf1.prof[,5]),lty=2,col=2)
lines(M, -tf1.prof[,6]-min(-tf1.prof[,6]),lty=5,col=2)
lines(M, -tf1.prof[,7]-min(-tf1.prof[,7]),lty=1,col=4)
lines(M, tf1.prof[,8]-min(tf1.prof[,8]),lty=2,col=4)
lines(M, -tf1.prof[,9]-min(-tf1.prof[,9]),lty=5,col=4)
legend("topleft", bty ='n', lwd=1, lty=c(1,5,1,5,1,2,5),col=c(1,1,2,2,4,4,4),
legend = c('Trawl survey','Commercial cpue','Trawl size','Winter Pot size','Summer Commercial size','recruits','Observer size') )







		