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
ny <- dat$ny
t1 <- dat$ny-1
dat$ny <- dat$ny-1
dat$sc <- dat$sc[-ny]
dat$tc <- dat$tc[-ny]
dat$te <- dat$te[-ny]
dat$stcpue <- dat$stcpue[-ny]
dat$stcpuese <- dat$stcpuese[-ny]
dat$ys <- dat$ys[-ny]
dat$onc <- dat$onc[,-ny]
dat$ooc <- dat$ooc[,-ny]
dat$twc <- dat$twc[-ny]
dat$tws <- dat$tws[-ny]
dat$twst <- dat$twst[-ny]
dat$tsc <- dat$tsc[-ny]
dat$tsct <- dat$tsct[-ny]
if (dat$it[length(dat$it)]==ny){
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
if (dat$iw[length(dat$iw)]==ny){
	y <- length(dat$iw)
	dat$nyw <- dat$nyw -1
	dat$sw <- dat$sw[-y]
    dat$wcpue <- dat$wcpue[-y]
	dat$iw <- dat$iw[-y]
	dat$onw <- dat$onw[,-y]
	dat$oow <- dat$oow[,-y]
    }
if (dat$io[length(dat$io)]==ny){
	z <- length(dat$io)
	dat$nyo <- dat$nyo -1
	dat$io <- dat$io[-z]
    dat$so <- dat$so[-z]
	dat$ono <- dat$ono[,-z]
	dat$ooo <- dat$ooo[,-z]
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
