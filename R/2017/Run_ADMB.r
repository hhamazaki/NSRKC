############################################################################
#  Run ADMB from Anywhere
############################################################################
rm(list=ls(all=TRUE))
library(PBSadmb)

############################################################################
#   1.0 Define ADMB exe file directry and name  
############################################################################
admodeldir <- paste('C:/Users/sid/Documents/MayCPT/')
admodel <- paste('DutchSc1')
setwd(admodeldir)

############################################################################
#   1.1 Define data file directry and name  
############################################################################
# infun: input file name 
datadir <- paste('C:/Users/sid/Documents/MayCPT/')
dataname <- paste('DutchSc1.dat')

############################################################################
#   1.2  Create model running environment  
############################################################################
# ADMB model file name
prefix <- paste0(admodeldir,admodel)
# out: ADMB parameter ouput file name 
out.rep <- paste0(prefix,'.rep')
out.std <- paste0(prefix,'.std')
out.par <- paste0(prefix,'.std')
# ADMB data model file name
admbdata <- paste0(datadir,dataname)
# Define ADMB argument 
argvec <- paste(' -ind',admbdata)
runAD(prefix=prefix, argvec = argvec, logfile=TRUE, add = FALSE, verbose=TRUE)


# Read output file
#rept <- reptoRlist(out.rep)
#rept
#par1 <- read.table(out.std, sep="",header=T)
#pars <- data.frame(t(par1[,3]))
#names(pars) <- t(par1[,2])

