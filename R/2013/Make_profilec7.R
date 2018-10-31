###############################################################################################
#   Crab Data profile likelihood construction for 2013 
#   This program read original data, modify save, run ADMB get likelihood value, and repeat this.
###############################################################################################
library(PBSadmb)
setwd('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2013/Model')

###############################################################################################
#   Crab Data construction 
###############################################################################################
# Crab data reading routine 
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/read_crab_data2013.R')
# Crab ouptutdata reading routine 
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/read_rep_data2013.R')
# set R-data file name 
fn <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/NSRKC2012R.txt')
# Read crab data into R List file 
dat <- datatoR(fn)

###############################################################################################
#   set output 
###############################################################################################
nrep <- 20
M.prof <- rep(0,nrep)
tf1.prof <- matrix(0,ncol=15,nrow=nrep)
s1 <- 0
incl <- 10

for (i in 1:nrep){
#dat$ms <- 1
#dat$M <- 0.22
#dat$slt6 <- 1
dat$maxss <- s1 + incl*i
file.remove(paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/','ns3n2012R_M.dat',sep=""))
out_file <- file(paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/','ns3n2012R_M.dat',sep=""), open='a')
for (j in seq_along(dat)){
    write.table(paste("#",names(dat)[j]), file=out_file, sep="", dec=".", quote=FALSE, col.names=FALSE, row.names=FALSE)  
	write(t(dat[[j]]), file=out_file, sep=" ", 
	ncolumns= ifelse(is.numeric(dim(dat[[j]])[2]),dim(dat[[j]])[2],length(dat[[j]])))
	}
close(out_file) 

# Run the ADMB 
argvec <- paste(' -ind ','ns3n2012R_M.dat',sep="")
runAD(prefix="ns3n2011c7", argvec = argvec, verbose=FALSE)
fn <- 'C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/ns3n2011c7.rep' 
# Run data read routin
rept <- reptoRlist(fn)
M.prof[i] <- as.double(scan(fn,skip=108,nlines=1,quiet=T,what="")) 
tf1.prof[i,] <- as.double(scan(fn,skip=110,nlines=1,quiet=T,what="")) 
}
M.prof
tf1.prof


# Clean up data 
#M.prof <- ifelse(M.prof>12000,NA,M.prof)
#tf1.prof[,2] <- ifelse(tf1.prof[,2]>40,NA,tf1.prof[,2])
#tf1.prof[,3] <- ifelse(tf1.prof[,3]>7,NA,tf1.prof[,3])
#tf1.prof[,5] <- ifelse(tf1.prof[,5]>6,NA,tf1.prof[,5])
#tf1.prof[,6] <- ifelse(-tf1.prof[,6]>2700,NA,tf1.prof[,6])
#tf1.prof[,7] <- ifelse(-tf1.prof[,7]>1300,NA,tf1.prof[,7])
#tf1.prof[,9] <- ifelse(-tf1.prof[,9]>2850,NA,tf1.prof[,9])
#tf1.prof[,11] <- ifelse(-tf1.prof[,11]>3700,NA,tf1.prof[,11])
#tf1.prof[,13] <- ifelse((tf1.prof[,13]>0.7|tf1.prof[,13]<0.2),NA,tf1.prof[,13])
#tf1.prof[,15] <- ifelse(-tf1.prof[,15]>550,NA,tf1.prof[,15])


par(mfrow=c(4,3),mar = c(2, 2, 2, 2))
M <- seq(s1+incl, s1+nrep*incl, by=incl)
plot(M, M.prof,type='l', main='Total negative log likelihood')
plot(M, tf1.prof[,2],type='l', main='Trawl survey')
plot(M, tf1.prof[,3],type='l', main='Pot survey')
plot(M, tf1.prof[,5],type='l', main='Commercial effort')
plot(M, -tf1.prof[,6],type='l',main='Trawl size')
plot(M, -tf1.prof[,7],type='l',main='Summer Pot size')
plot(M, -tf1.prof[,9],type='l',main='Winter Pot size')
plot(M, -tf1.prof[,11],type='l',main='Summer Commercial size')
plot(M, tf1.prof[,13],type='l', main='recruits')
plot(M, -tf1.prof[,15],type='l',main='Observer size')


par(mfrow=c(4,3),mar = c(2, 2, 2, 2))
M <- seq(s1+incl, s1+nrep*incl, by=incl)
plot(M, M.prof,type='l', ylim=c(10980,11050),main='Total negative log likelihood')
plot(M, tf1.prof[,2],type='l',ylim=c(24,40), main='Trawl survey')
plot(M, tf1.prof[,3],type='l', ylim=c(1.2,4.0),main='Pot survey')
plot(M, tf1.prof[,5],type='l', ylim=c(4.5,6.0),main='Commercial effort')
plot(M, -tf1.prof[,6],type='l',ylim=c(2580,2610),main='Trawl size')
plot(M, -tf1.prof[,7],type='l',ylim=c(1272,1290),main='Summer Pot size')
plot(M, -tf1.prof[,9],type='l',ylim=c(2820,2840),main='Winter Pot size')
plot(M, -tf1.prof[,11],type='l',ylim=c(3620,3700),main='Summer Commercial size')
plot(M, tf1.prof[,13],type='l',ylim=c(0.30,0.80), main='recruits')
plot(M, -tf1.prof[,15],type='l',ylim=c(545,548),main='Observer size')





		