###############################################################################################
#   Reduced Retrospective Crab Data construction for past 10 years  
###############################################################################################
# Reduce data year by year and save as different file. 
makeretrodata <- function(dat,rn,outfn){
# Save original file as filename_0.dat
file.remove(paste(outfn,0,'.dat',sep=""))
makeADdata(dat,paste(outfn,0,'.dat',sep=""))
 # Create reduced data file and save as filenme_i.dat
for (i in 1:rn){
ny <- dat$lyear - dat$fyear + 1
dat$lyear <- dat$lyear -1
dat$tc <- dat$tc[-ny]
dat$te <- dat$te[-ny]
dat$stcpue <- dat$stcpue[-ny]
dat$secpue <- dat$secpue[-ny]
dat$ys <- dat$ys[-ny]
dat$nonc <- dat$nonc[,-ny]
dat$nooc <- dat$nooc[,-ny]
dat$twc <- dat$twc[-ny]
dat$tws <- dat$tws[-ny]
dat$twst <- dat$twst[-ny]
dat$tsc <- dat$tsc[-ny]
dat$tsct <- dat$tsct[-ny]
# Trawl Survey
if (dat$it[length(dat$it)]==(dat$lyear+1)){
	x <- length(dat$it[dat$it < dat$lyear+1])
	n <- length(dat$it) - x
	dat$nyt <- dat$nyt -n
	dat$it <- dat$it[1:x]
	dat$tt <- dat$tt[1:x]
	dat$pct <- dat$pct[1:x]
	dat$yt <- dat$yt[1:x]
	dat$cv <- dat$cv[1:x]
	dat$nont <- dat$nont[,1:x]
	dat$noot <- dat$noot[,1:x]	
}
# Winter Pot Survey
if (dat$iw[length(dat$iw)]==(dat$lyear+1)){
	y <- length(dat$iw)
	dat$nyw <- dat$nyw -1
    dat$wcpue <- dat$wcpue[-y]
	dat$iw <- dat$iw[-y]
	dat$nonw <- dat$nonw[,-y]
	dat$noow <- dat$noow[,-y]
    }
# Observer Discards 
if (dat$iod[length(dat$iod)]==(dat$lyear+1)){
	z <- length(dat$iod)
	dat$nyod <- dat$nyod -1
	dat$iod <- dat$iod[-z]
	dat$nonod <- dat$nonod[,-z]
	dat$noood <- dat$noood[,-z]
    }
# Observer Total
if (dat$iot[length(dat$iot)]==(dat$lyear+1)){
	z <- length(dat$iot)
	dat$nyot <- dat$nyot -1
	dat$iot <- dat$iot[-z]
	dat$nonot <- dat$nonot[,-z]
	dat$nooot <- dat$nooot[,-z]
    }
# Observer Winter
if (dat$icw[length(dat$icw)]==(dat$lyear+1)){
	z <- length(dat$icw)
	dat$nywc <- dat$nywc -1
	dat$icw <- dat$icw[-z]
	dat$noncw <- dat$noncw[,-z]
	dat$noocw <- dat$noocw[,-z]
    }	
	
file.remove(paste(outfn,i,'.dat',sep=""))
makeADdata(dat,paste(outfn,i,'.dat',sep=""))
}
}
