library(R2admb)
library(ggplot2)
library(PBSadmb)
###############################################################################################
#   Crab Data construction 
###############################################################################################
# Read reproRlist.R for program use
# Crab data reading routine 
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/read_crab_data2013.R')
# Read crab data into R List file 
dat <- datatoR('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2013/Model/ns3n2012R.txt')
# Read Make Crab retro data reading routine 
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/make_retro_data2013.R')
# Make Crab retro data reading routine 
outfn <-'C:/Projects/Norton_Sound/NSCrab/Prediction_model/2013/Test1/ns3n2013R_'
makeretrodata(dat,2,outfn)
# Read Make Crab forward data reading routine 
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/make_forward_data2013.R')
# Make Crab retro data reading routine 
outfn <-'C:/Projects/Norton_Sound/NSCrab/Prediction_model/2013/Test2/ns3n2013R_'
makeforwarddata(dat,2,outfn)

###############################################################################################
#   Run ADMB model and Extract output to R data 
###############################################################################################

# Create empty data 
MMB <- matrix(0,ncol=11, nrow=37)
# Run the ADMB 
for (i in 0:10){
argvec <- paste(' -ind',outfn,i,'.dat',sep="")
runAD(prefix="ns3n31", argvec = argvec,add=False)
fn <- 'C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/ns3n31.rep'
MMB[,i+1] <- as.double(scan(fn,skip=85,nlines=1,quiet=T,what="",nmax=37)) 
}
MMB

year <- seq(1976,2012)
plot(year, MMB[,1], type ='l')
for (i in 1:10){
a <- 37-i
lines(year[1:a],MMB[1:a,i+1], col = i)
}


		