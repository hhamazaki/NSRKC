###############################################################################################
#   Read ADMB dat Data to R
###############################################################################################
	
datatoR = function(fn,par=F)
{
if (par) flush=F else flush=T
ifile=scan(fn,what="character",flush=flush, blank.lines.skip=F,quiet=T)

#if(par) ifile=ifile[-c(2:3,5,7:8,10,12,14:15)][ifile[-c(2:3,5,7:8,10,12,14:15)]!='#']
idx=sapply(as.double(ifile),is.na)
vnam=ifile[idx] #list names
nv=length(vnam) #number of objects
A=list()
ir=0
for(i in 1:nv)
{
	ir=match(vnam[i],ifile)
	if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
	dum=NA
	if(par) dum=as.numeric(as.character(ifile[(ir+1):(irr-1)])) else {
		if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=T,what="")) 
		if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=T))
		}

	if(is.numeric(dum))#Logical test to ensure dealing with numbers
	name=sub('#','',vnam[i])
	A[[name]]=dum
}
return(A)
}


###############################################################################################
#   Make R data to ADMB Data
###############################################################################################
	
makeADdata <- function(dat,outfn){
# Save original file as filename_0.dat
file.remove(paste(outfn))
out_file <- file(paste(outfn), open='a')
for (j in seq_along(dat)){
    write.table(paste0("#",names(dat)[j]), file=out_file, sep="", dec=".", quote=FALSE, col.names=FALSE, row.names=FALSE)  
	write(t(dat[[j]]), file=out_file, sep=" ", 
	ncolumns= ifelse(is.numeric(dim(dat[[j]])[2]),dim(dat[[j]])[2],length(dat[[j]])))
	}
close(out_file) 
}

###############################################################################################
#   Reac Crab report Data to R
###############################################################################################
# Define the report data
reptoRlist = function(fn,par=F)
{
if (par) flush=F else flush=T
ifile=scan(fn,what="character",flush=flush, blank.lines.skip=F,quiet=T)
idx=sapply(as.double(ifile),is.na)
vnam=ifile[idx] #list names
nv=length(vnam) #number of objects
A=list()
ir=0
for(i in 1:nv)
{
	ir=match(vnam[i],ifile)
	if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
	dum=NA
	if(par) dum=as.numeric(as.character(ifile[(ir+1):(irr-1)])) else {
		if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=T,what="")) 
		if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=T))
		}

	if(is.numeric(dum))#Logical test to ensure dealing with numbers
	name=sub(':','',vnam[i])
	A[[name]]=dum
}
return(A)
}

