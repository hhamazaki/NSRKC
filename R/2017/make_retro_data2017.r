###############################################################################################
#   Reduced Retrospective Crab Data construction for past 10 years  
###############################################################################################
# Reduce data year by year and save as different file. 
makeretrodata <- function(dat,rn,outfn){
# Save original file as filename_0.dat
file.remove(paste(outfn,0,'.dat',sep=""))
out_file <- file(paste(outfn,0,'.dat',sep=""), open='a')
for (j in seq_along(dat)){
    write.table(paste("#",names(dat)[j]), file=out_file, sep="", dec=".", quote=FALSE, col.names=FALSE, row.names=FALSE)  
	write(t(dat[[j]]), file=out_file, sep=" ", 
	ncolumns= ifelse(is.numeric(dim(dat[[j]])[2]),dim(dat[[j]])[2],length(dat[[j]])))
	}
close(out_file) 
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
if (dat$it[length(dat$it)]==(dat$lyear+1)){
	x <- length(dat$it)
	dat$nyt <- dat$nyt -1
	dat$it <- dat$it[-x]
	dat$tt <- dat$tt[-x]
	dat$pct <- dat$pct[-x]
	dat$yt <- dat$yt[-x]
	dat$cv <- dat$cv[-x]
	dat$nont <- dat$nont[,-x]
	dat$noot <- dat$noot[,-x]	
}
if (dat$iw[length(dat$iw)]==(dat$lyear+1)){
	y <- length(dat$iw)
	dat$nyw <- dat$nyw -1
    dat$wcpue <- dat$wcpue[-y]
	dat$iw <- dat$iw[-y]
	dat$nonw <- dat$nonw[,-y]
	dat$noow <- dat$noow[,-y]
    }
if (dat$io[length(dat$io)]==(dat$lyear+1)){
	z <- length(dat$io)
	dat$nyo <- dat$nyo -1
	dat$io <- dat$io[-z]
	dat$nono <- dat$nono[,-z]
	dat$nooo <- dat$nooo[,-z]
    }

file.remove(paste(outfn,i,'.dat',sep=""))
out_file <- file(paste(outfn,i,'.dat',sep=""), open='a')
for (j in seq_along(dat)){
    write.table(paste("#",names(dat)[j]), file=out_file, sep="", dec=".", quote=FALSE, col.names=FALSE, row.names=FALSE)  
	write(t(dat[[j]]), file=out_file, sep=" ", 
	ncolumns= ifelse(is.numeric(dim(dat[[j]])[2]),dim(dat[[j]])[2],length(dat[[j]])))
	}
close(out_file) 
}
}
