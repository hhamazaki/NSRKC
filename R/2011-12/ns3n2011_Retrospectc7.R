library(PBSadmb)
library(gplots)
setwd('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model')
###############################################################################################
#   Crab Data construction 
###############################################################################################

###############################################################################################
#   1.0 Read source file 
###############################################################################################
# Read read.R for program use
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/read_crab_data.R')
# Make retro data
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/make_retro_data.R')

###############################################################################################
#   1.1 Define inpug file name and output file name and locations. 
###############################################################################################
# infun: input file name 
infn <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/NSRKC2012R.txt')
# outfn: output file name  
outfn <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Test1/ns3n2011R')

###############################################################################################
#   1.2 Red ADMB data into R  
###############################################################################################
# Read crab data into R List file 
dat <- datatoR(infn)

# insert sta tements to change data
# dat$M <- 0.32
 dat$ms <- 3.6
# dat$slt6 <- 1.0
 dat$lamc <- 100
 dat$maxss <- 50
# dat$efn <-0.5
lt <- 'lamc=100, maxss = 50 ms6 =3.6'
###############################################################################################
#   1.3 Create retrospective data  
###############################################################################################
# Determine the number of years to make retrospect
rn <- 15
ny <- dat$ny  # the number of years modeled
# The file will be written from 0 to 10
makeretrodata(dat,rn,outfn)

###############################################################################################
#   1.4  Run ADMB model and Extract output to R data 
###############################################################################################
# Make sure change data 
# Create empty data 
legal <- matrix(0,ncol=rn+1, nrow=ny)
f <- numeric(rn+1)
# set moedel prefix
prefix <- paste('ns3n2011c7')
# set retro data prefix
datfn <- paste('ns3n2011R_')
# Run the ADMB 
for (i in 0:rn){
argvec <- paste(' -ind ',outfn,i,'.dat',sep="")
runAD(prefix=prefix, argvec = argvec, logfile=FALSE, add = FALSE, verbose=FALSE)
# Read data from the output 
fn <- 'C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/ns3n2011c7.rep'
legal[,i+1] <- as.double(scan(fn,skip=94,nlines=1,quiet=T,what="",nmax=ny)) 
f[i+1] <- as.double(scan(fn,skip=108,nlines=1,quiet=T,what="",nmax=ny)) 
fnp <- 'C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/ns3n2011c7.std'
#file.rename(paste(fn), paste(outfn,i,'.rep',sep=""))
#file.rename(paste(fnp), paste(outfn,i,'.std',sep=""))
}
legal
f

###############################################################################################
#   1.5  plot 
###############################################################################################
windows(record=TRUE)
year <- seq(1976,1975+ny)
plot(year, legal[,1], type ='l')
for (i in 1:rn){
a <- ny-i
lines(year[1:a],legal[1:a,i+1], col = i)
}

# Extract hind cast data for past 10 years
hind <- legal[(ny-rn):(ny-1),1]
hindpred <- numeric(rn)
for (i in 1:rn){
hindpred[rn+1-i] <- legal[ny-i,1+i]
}
hindpred
minx <- min(min(hindpred),min(hind))
maxx <- max(max(hindpred),max(hind))
year2 <- seq(1975+ny-rn,1976+ny-2)
plot(year2, hind, type = 'l', ylim=c(minx,maxx), ylab = 'Legal abundance', xlab = 'Year')
lines(year2, hindpred, col = 2)
# Estimate legal abundance
plegal <- colSums((dat$ont+dat$oot)*dat$lg)
oblegal <- dat$tt*plegal
points(year[dat$it],oblegal)
ciw <- 2*oblegal*dat$cv3
plotCI(year[dat$it],oblegal,uiw=ciw,ylab = 'Legal abundance', xlab = 'Year',ylim=c(0,6000))
lines(year,legal[,1])
lines(year2,hindpred,col=2)

plotCI(year[dat$it],oblegal,uiw=ciw,ylab = 'Legal abundance', xlab = 'Year',ylim=c(0,3000),xlim=c(1975+ny-rn,1976+ny-2))
lines(year,legal[,1],lwd=2)
lines(year2,hindpred,col=2,lwd=2)
legend("topright",lty =c(0,1,1,0), lwd =c(0,2,2),col=c(NA,2,1,1), pch=c(NA,NA,NA,1), bty='n', legend=c(lt,'Hindcast predicted','Hindcast 2011','Observed'))

plotCI(year[dat$it],oblegal,uiw=ciw,ylab = 'Legal abundance', xlab = 'Year',ylim=c(0,6000))
for (i in 1:rn){
a <- ny-i
lines(year[1:a],legal[1:a,i+1], col = i)
}
lines(year,legal[,1],lwd=2)
lines(year2,hindpred,col=2,lwd=2)
legend("topright",lty=c(0,1,1,0),lwd =c(0,2,2), col=c(NA,2,1,1), pch=c(NA,NA,NA,1), bty='n', legend=c(lt,'Hindcast predicted','Hindcast 2011','Observed'))


# Calcurate retrospective bias
phind <- (hindpred-hind)/hind
avhind <- mean(phind)
avhind
coef(lm(hind ~ -1 + hindpred))




