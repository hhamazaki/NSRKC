###############################################################################################
#   Forward Retrospective Crab Data construction for past 10 years  
###############################################################################################
# Reduce data year by year and save as different file. 
makeforwarddata <- function(dat,rn,outfn){
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
# Change commercial fishery years
dat$scy1 <- dat$scy1 -1
dat$scy2 <- dat$scy2 -1
dat$scp <- dat$scp -1
dat$ny <- dat$ny-1
# Remove data
dat$sc <- dat$sc[-1]
dat$tc <- dat$tc[-1]
dat$te <- dat$te[-1]
dat$stcpue <- dat$stcpue[-1]
dat$stcpuese <- dat$stcpuese[-1]
dat$ys <- dat$ys[-1]
dat$onc <- dat$onc[,-1]
dat$ooc <- dat$ooc[,-1]
dat$twc <- dat$twc[-1]
dat$tsc <- dat$tsc[-1]

# Change Trawl Data
dat$it <- dat$it - 1
if (dat$it[1] == 0){
	dat$it <- dat$it[-1]
	dat$nyt <- dat$nyt-1
	dat$tt <- dat$tt[-1]
	dat$pct <- dat$pct[-1]
	dat$yt <- dat$yt[-1]
	dat$st <- dat$st[-1]
	dat$cv1 <- dat$cv1[-1]
	dat$cv2 <- dat$cv2[-1]
	dat$cv3 <- dat$cv3[-1]
	dat$ont <- dat$ont[,-1]
	dat$oot <- dat$oot[,-1]	
}

# Change Winter Survey Data
dat$iw <- dat$iw - 1
if (dat$iw[1] == 0){
	dat$nyw <- dat$nyw-1
	dat$sw <- dat$sw[-1]
	dat$iw <- dat$iw[-1]
	dat$onw <- dat$onw[,-1]
	dat$oow <- dat$oow[,-1]
}
dat$ipre <- dat$ipre -1
dat$io <- dat$io -1

# Change Pot Survey Data
dat$ip <- dat$ip -1
if (dat$ip[1] == 0){
	dat$nyp <- dat$nyp -1
	dat$ip <- dat$ip[-1]
	dat$tp <- dat$tp[-1]
	dat$yp <- dat$yp[-1]
	dat$sp <- dat$sp[-1]
	dat$onp <- dat$onp[,-1]
	dat$oop <- dat$oop[,-1]
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
