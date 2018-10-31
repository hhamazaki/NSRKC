library(R2admb)
library(ggplot2)
library(PBSadmb)
###############################################################################################
#   Crab Data construction 
###############################################################################################

# Read reproRlist.R for program use
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/reptoRlist1.R')
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/Read_rep_data.r')
# Read crab data into R List file 
dat=reptoRlist('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/ns3n2011R.txt')
# Change M to 0.32
# dat$M <- 0.32
# Output file 
file.remove(paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/','ns3n2011R_',0,'.dat',sep=""))
out_file <- file(paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/','ns3n2011R_',0,'.dat',sep=""), open='a')
for (j in seq_along(dat)){
    write.table(paste("#",names(dat)[j]), file=out_file, sep="", dec=".", quote=FALSE, col.names=FALSE, row.names=FALSE)  
	write(t(dat[[j]]), file=out_file, sep=" ", 
	ncolumns= ifelse(is.numeric(dim(dat[[j]])[2]),dim(dat[[j]])[2],length(dat[[j]])))
	}
close(out_file) 
fn <- 'C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/ns3n31c2.rep'
rept <- readReport(fn)

###############################################################################################
#   Reduced Retrospective Crab Data construction for past 10 years  
###############################################################################################
# Reduce data year by year and save as different file. 
for (i in 1:10){
t1 <- dat$ny-1
dat$ny <- dat$ny-1
dat$sc <- dat$sc[-t1]
dat$tc <- dat$tc[-t1]
dat$te <- dat$te[-t1]
dat$ys <- dat$ys[-t1]
dat$onc <- dat$onc[,-t1]
dat$ooc <- dat$ooc[,-t1]
dat$twc <- dat$twc[-t1]
dat$tsc <- dat$tsc[-t1]
if (dat$it[length(dat$it)]==t1){
	x <- length(dat$it)
	dat$nyt <- dat$nyt-1
	dat$it <- dat$it[-x]
	dat$tt <- dat$tt[-x]
	dat$pct <- dat$pct[-x]
	dat$yt <- dat$yt[-x]
	dat$st <- dat$st[-x]
	dat$cv1 <- dat$cv1[-x]
	dat$cv2 <- dat$cv2[-x]
	dat$cv3 <- dat$cv3[-x]
	dat$ont <- dat$ont[,-x]
	dat$oot <- dat$oot[,-x]	
}
if (dat$iw[length(dat$iw)]==t1){
	y <- length(dat$iw)
	dat$nyw <- dat$nyw -1
	dat$sw <- dat$sw[-y]
	dat$iw <- dat$iw[-y]
	dat$onw <- dat$onw[,-y]
	dat$oow <- dat$oow[,-y]
}
file.remove(paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/','ns3n2011R_',i,'.dat',sep=""))
out_file <- file(paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/','ns3n2011R_',i,'.dat',sep=""), open='a')
for (j in seq_along(dat)){
    write.table(paste("#",names(dat)[j]), file=out_file, sep="", dec=".", quote=FALSE, col.names=FALSE, row.names=FALSE)  
	write(t(dat[[j]]), file=out_file, sep=" ", 
	ncolumns= ifelse(is.numeric(dim(dat[[j]])[2]),dim(dat[[j]])[2],length(dat[[j]])))
	}
close(out_file) 
}

###############################################################################################
#   Run ADMB model and Extract output to R data 
###############################################################################################

# Create empty data 
MMB <- matrix(0,ncol=11, nrow=37)
# Run the ADMB 
for (i in 0:10){
argvec <- paste(' -ind ','ns3n2011R_',i,'.dat',sep="")
runAD(prefix="ns3n31", argvec = argvec)
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


		