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
outfn <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Test1/ns3n2011model.dat')

###############################################################################################
#   1.2 Red ADMB data into R  
###############################################################################################
# Read crab data into R List file 
dat <- datatoR(infn)

makedata <- function(dat,outfn){
# Save original file as filename_0.dat
file.remove(paste(outfn))
out_file <- file(paste(outfn), open='a')
for (j in seq_along(dat)){
    write.table(paste("#",names(dat)[j]), file=out_file, sep="", dec=".", quote=FALSE, col.names=FALSE, row.names=FALSE)  
	write(t(dat[[j]]), file=out_file, sep=" ", 
	ncolumns= ifelse(is.numeric(dim(dat[[j]])[2]),dim(dat[[j]])[2],length(dat[[j]])))
	}
close(out_file) 
}
lt <- c(
	"Baseline 2012",
	"lamc = 100: Model 4",
	"maxss = 10: Model 7",
	"M = 0.32: Model 8",
	"M = 0.4 ms6 = 1.0: Model 10",
	"maxss = 50, lamc = 50: Model 11",
	"maxss = 50, lamc = 50, ms6 = 3.6: Model 12"
	)

tasks <- function(dat, itasks){
	if(itasks == 0) makedata(dat,outfn)
	if(itasks == 1) {
	dat <- datatoR(infn)
	dat$lamc <- 100
	makedata(dat,outfn)
	}
	if(itasks == 2) {
	dat <- datatoR(infn)
	dat$maxss <- 10
	makedata(dat,outfn)
	}
	if(itasks == 3) {
	dat <- datatoR(infn)
	dat$M <- 0.32
	makedata(dat,outfn)
	}
	if(itasks == 4) {
	dat <- datatoR(infn)
	dat$M = 0.4
	dat$ms6 <- 1.0
	makedata(dat,outfn)
	}
	if(itasks == 5) {
	dat <- datatoR(infn)
	dat$lamc <- 50
	dat$maxss <- 50
	makedata(dat,outfn)
	}
	if(itasks == 6) {
	dat <- datatoR(infn)
	dat$lamc <- 50
	dat$maxss <- 50
	dat$ms <-3.6
	makedata(dat,outfn)
	}
	}
###############################################################################################
#   1.3 Create retrospective data  
###############################################################################################
# Determine the number of years to make retrospect
rn <- 6
ny <- dat$ny
###############################################################################################
#   1.4  Run ADMB model and Extract output to R data 
###############################################################################################
# Make sure change data 
# Create empty data 
legal <- matrix(0,ncol=rn+1, nrow=ny)
mmb <- matrix(0,ncol=rn+1, nrow=ny)
f <- numeric(rn+1)

# set moedel prefix
prefix <- paste('ns3n2011c7')
# set retro data prefix
datfn <- paste('ns3n2011M_')
# Run the ADMB 
for (i in 0:rn){
tasks(dat, i)
argvec <- paste(' -ind ',outfn,sep="")
runAD(prefix=prefix, argvec = argvec, logfile=FALSE, add = FALSE, verbose=FALSE)
# Read data from the output 
fn <- 'C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/ns3n2011c7.rep'
legal[,i+1] <- as.double(scan(fn,skip=94,nlines=1,quiet=T,what="",nmax=ny)) 
mmb[,i+1] <- as.double(scan(fn,skip=96,nlines=1,quiet=T,what="",nmax=ny)) 
f[i+1] <- as.double(scan(fn,skip=108,nlines=1,quiet=T,what="",nmax=ny)) 

fnp <- 'C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/ns3n2011c7.std'
}
legal
mmb
f

###############################################################################################
#   1.5  plot 
###############################################################################################

windows(record=TRUE)
year <- seq(1976,1975+ny)
# Estimate legal abundance
plegal <- colSums((dat$ont+dat$oot)*dat$lg)
oblegal <- dat$tt*plegal
ciw <- 2*oblegal*dat$cv3
plotCI(year[dat$it],oblegal,uiw=ciw,ylab = 'Legal abundance', xlab = 'Year',ylim=c(0,7000))
for (i in 1:rn){
lines(year,legal[,i+1], col = ifelse(i+1<6,1,2), lty = i+1)
}
lines(year,legal[,1],lwd=2)
legend(1985, 6000, bty ='n', lwd=2, col = 1, "Baseline 2012")
for (i in 1:rn){
legend(1985, 6000-200*i, bty ='n', lwd=1, lty=i+1, col = ifelse(i+1 < 6,1,2), lt[i+1])
}
plot(year,mmb[,1],lwd=2,type = 'l', ylab = 'MMB', xlab = 'Year',ylim=c(0,20000))
for (i in 1:rn){
lines(year,mmb[,i+1], col = ifelse(i+1<6,1,2), lty = i+1)
}
legend(1985, 12000, bty ='n', lwd=2, col = 1, "Baseline 2012")
for (i in 1:rn){
legend(1985, 12000-600*i, bty ='n', lwd=1, lty=i+1, col = ifelse(i+1 < 6,1,2), lt[i+1])
}


