##############################################################################################################
#  Figure producing R function rutine
##############################################################################################################
library(gplots)
library(PBSadmb)
setwd('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model')
# Read reproRlist.R for program use
# Crab data reading routine 
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/read_crab_data.R')
# Crab ouptutdata reading routine 
source('C:/Projects/Norton_Sound/NSCrab/Prediction_model/Models/R/read_rep_data2.R')
# Read crab data into R List file 
dat <- datatoR('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/NSRKC2012R.txt')
# outfn: output file name  
outfn <- paste('C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Test1/ns3n2011model.dat')
# Determine the output file location
fn <- 'C:/Projects/Norton_Sound/NSCrab/Prediction_model/2012/Model/ns3n2011c7.rep'
prefix <- paste('ns3n2011c7')
###############################################################################################
#   1.2 Red ADMB data into R  
###############################################################################################
# Read crab data into R List file 
# insert sta tements to change data
# dat$M <- 0.32
  dat$ms <- 3.6
 dat$lamc <- 50
 dat$maxss <- 50
# dat$efn <-0.5

######### Remove the last winter pot survey data ##############################################
#	y <- length(dat$iw)
#	dat$nyw <- dat$nyw -1
#	dat$sw <- dat$sw[-y]
#	dat$iw <- dat$iw[-y]
#	dat$onw <- dat$onw[,-y]
#	dat$oow <- dat$oow[,-y]
###############################################################################################	
	
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
makedata(dat,outfn)

# set moedel prefix

# Run the ADMB 
argvec <- paste(' -ind ',outfn,sep="")
runAD(prefix=prefix, argvec = argvec, logfile=FALSE, add = FALSE, verbose=FALSE)


# Read output file
rept <- reptoRlist(fn)
rept

##############################################################################################################
#  Calculate Effective sample size
##############################################################################################################
#trawl survey effective sample size
ot <- dat$ont+dat$oot
tefn <-colSums(rept$ent*(1-rept$ent))/colSums((ot-rept$ent)^2)
tefn

#Summer pot survey effective sample size
spp <- rept$enp+rept$eop
sp <- dat$onp+dat$oop
spefn <-colSums(spp*(1-spp))/colSums((spp-sp)^2)

#commercial fishery effective sample size
scp <- rept$enc+rept$eoc
sc <- dat$onc+dat$ooc
scefn <-colSums(scp*(1-scp))/colSums((scp-sc)^2)
scefn <-scefn[-c(1,16)]
#Winter Pot survey effective sample size
wpp <- rept$enw+rept$eow
wp <- dat$onw+dat$oow
wpefn <-colSums(wpp*(1-wpp))/colSums((wpp-wp)^2)

# Observer survey effective sample size
opp <- rept$eno+rept$eoo
op <- dat$ono+dat$ooo
opefn <-colSums(opp*(1-opp))/colSums((opp-op)^2)
tefn
spefn
scefn
wpefn
opefn






# Create year data 
year <- seq(1976,1975+dat$ny)
lenclass <- c('74-83','84-93','94-103','104-113','114-123','>123')

windows(record=TRUE)
##############################################################################################################
#  Figure for trends.
##############################################################################################################

# Figure 1. Total estimated abundance 
total.abundance <- colSums(rept$nps+rept$ops)/1000
maxlim <- max(total.abundance)
# Plot trawl survey data with CI
ciw <- 2*dat$tt*dat$cv3/1000
plotCI(year[dat$it],dat$tt/1000,uiw=ciw,ylab = 'Total Crab Abundance (million)', xlab = 'Year', ylim= c(0,maxlim),bty='l')
lines(year,total.abundance[1:dat$ny],lwd=2)
legend("topright", lty = c(1,0), lwd = c(2,0), pch = c(NA,1), legend = c('Predicted','Observed'))
title(main='Trawl survey crab abundance')
sigmat <- sqrt(sum((log(dat$tt)-log(rept$ett1+rept$ett2))^2)/length(dat$tt))
sigmap <- sqrt(sum((log(dat$tp)-log(rept$etp1+rept$etp2))^2)/length(dat$tp))


# Figure 2. Total, Recruites, Legals estimated abundance 
plot(year, total.abundance[1:dat$ny], ylab = 'Abundance (million crabs)', typ = 'l', lwd=2, xlab = 'Year',ylim= c(0,maxlim),bty='l')
lines(year, rept$Legal/1000, lty=2)
lines(year, rept$Rec/1000, lty=3)
legend("topright",lty = c(1,2,3), legend = c('total','legal','recruits'))
title(main='Modeled crab abundance')
# Figure 3. Effort data 
plot(year[-dat$ny], rept$ete[-dat$ny]/1000, ylab = 'Effort (1000 pot lifts)', lwd=2, typ = 'l', xlab = 'Year', ylim=c(0,40), bty='l')
points(year[-dat$ny], dat$te[-dat$ny]/1000, pch=19)
legend("topright",lty =c(1,0), lwd = c(2,0), pch=c(NA,19), legend = c('Predicted','Observed'))
title(main='Summer commercial catch effort')
sigmae <- sqrt(sum((log(dat$te[-c(1,16,dat$ny)])-log(rept$ete[-c(1,16,dat$ny)]))^2)/length(dat$te[-c(1,16,dat$ny)]))

# Figure 4. Catch data 
par(mar=c(5,4,4,5))
tcatch <- dat$tc[-dat$ny]+dat$twc+dat$tsc
hrate <- tcatch/rept$Legal[-dat$ny]
plot(year[-dat$ny], tcatch[-dat$ny]/1000, ylab = 'Total Catch (million)', lwd=2, typ = 'l', xlab = 'Year', bty='u')
par(new=TRUE)
plot(year[-dat$ny],hrate,lty=2, lwd=2, typ = 'l',bty='l', xaxt='n', yaxt='n', xlab='',ylab='', ylim = c(0,0.5))
axis(4)
mtext("Predicted harvest rate",side=4,line=3)
legend("topright", lty = c(1,2), lwd = 2, legend=c('Total Catch', 'Estimated Harvest Rate'))
title(main='Total catch & Predicted harvest rate')



##############################################################################################################
#  Figure for size proportion 
##############################################################################################################
# Create commercial length vs predicted length freq figure 
par(mfrow=c(6,6),mar = c(1.2, 1.2, 1.2, 1.2),oma=c(3,1,2,1))
for (i in 2:15){
x <- plot(dat$onc[,i]+dat$ooc[,i],ylim=c(0,0.7),pch=19, cex=0.8, tck = -0.05, bty='l')
lines(rept$enc[,i]+rept$eoc[,i],lty=2)
text(2,0.6,paste(year[i]),cex=1.0)
}
for (i in 17:(dat$ny-1)){
x <- plot(dat$onc[,i]+dat$ooc[,i],ylim=c(0,0.7),pch=19, cex=0.8, tck = -0.05, bty='l')
lines(rept$enc[,i]+rept$eoc[,i],lty=2)
text(2,0.6,paste(year[i]),cex=1.0)
}
mtext("commercial harvest length: observed vs predicted", side=3, outer=T)
mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)

# Create winter length vs predicted length freq figure 
par(mfrow=c(6,6),mar = c(1.2, 1.2, 1.2, 1.2),oma=c(3,1,2,1))
for (i in 1:dat$nyw){
x <- plot(dat$onw[,i]+dat$oow[,i],ylim=c(0,0.5), pch=19, cex=0.8, tck = -0.05,bty='l', xlab = " ", ylab=" ")
lines(rept$enw[,i]+rept$eow[,i],lty=2)
text(2,0.45,paste(year[dat$iw[i]]),cex=1.0)
}
mtext("Winter pot length: observed vs predicted", side=3, outer=T)
mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)


# Create trawl length vs predicted length freq figure 
par(mfrow=c(6,6),mar = c(1.2, 1.2, 1.2, 1.2),oma=c(3,1,2,1))
for (i in 1:dat$nyt){
x <- plot(dat$ont[,i]+dat$oot[,i],ylim=(c(0,0.5)), pch=19,  cex=0.8, tck = -0.05, bty='l', xlab = " ", ylab=" ")
lines(rept$ent[,i],lty=2)
text(2,0.45,paste(year[dat$it[i]]),cex=1.0)
}
mtext("Trawl length: observed vs predicted", side=3, outer=T)
mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)
plot.new()
plot.new()
plot.new()
plot.new()
plot.new()

# Create pot survey length vs predicted length freq figure 
#par(mfrow=c(6,6),mar = c(1.2, 1.2, 1.2, 1.2),oma=c(3,1,2,1))
#par(mar = c(1.2, 1.2, 2, 1.2),oma=c(3,1,2,1))
for (i in 1:dat$nyp){
x <- plot(dat$onp[,i]+dat$oop[,i],ylim=(c(0,0.5)), pch=19,  cex=0.8, tck = -0.05, bty='l', xlab = " ", ylab=" ")
lines(rept$enp[,i]+rept$eop[,i],lty=2)
text(2,0.45,paste(year[dat$ip[i]]),cex=1.0)
}
mtext("Pot length: observed vs predicted", side=3, line = -25, outer=T)
#mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)
plot.new()
plot.new()
plot.new()
plot.new()
plot.new()
plot.new()
plot.new()
plot.new()
# Create observer survey length vs predicted length freq figure 
#par(mfrow=c(6,6),mar = c(1.2, 1.2, 1.2, 1.2),oma=c(3,1,2,1))
for (i in 1:dat$nyo){
x <- plot(dat$ono[,i]+dat$ooo[,i],ylim=(c(0,0.5)), pch=19,  cex=0.8, tck = -0.05, bty='l', xlab = " ", ylab=" ")
lines(rept$eno[,i]+rept$eoo[,i],lty=2)
text(2,0.45,paste(year[dat$io[i]]),cex=1.0)
}
mtext("Observer length: observed vs predicted", side=3, line = -40, outer=T)
#mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)


##############################################################################################################
#  Figure for bubble plot
##############################################################################################################
bubbleplot<-function(){
test2$radius <- 1.5*sqrt(abs(test2$fbdata)/pi) 
symbols(test2$years,test2$classes,circles=test2$radius, bg=ifelse(test2$fbdata<0,'white','black'),fg=ifelse(test2$fbdata==0,'white','black'), inches = FALSE, xlim=c(1976,1975+dat$ny-1),xlab = " ", ylab=" ")
}

par(mfrow=c(5,1), mar=c(2,2,2,1),oma=c(3,1,2,1))

# Create commercial length vs predicted length freq figure 
bdata <- dat$onc+dat$ooc - (rept$enc + rept$eoc)
fbdata <- as.vector(bdata)
classes <- rep(1:6,dat$ny)
years <- sort(rep(1976:(1975+dat$ny),6))
test <- data.frame(cbind(years,classes,fbdata))
test2 <- subset(test,((test$years != 1976)&(test$years != 1991)&(test$years != 2012)))
bubbleplot()
title(main='Commercial Harvest')

# Create Winter length vs predicted length freq figure 
bdata <- dat$onw+dat$oow - rept$enw - rept$eow
fbdata <- as.vector(bdata)
classes <- rep(1:6,dat$nyw)
years <-1
for (i in 1:dat$nyw){
years <- c(years,rep(year[dat$iw[i]],6))
}
years <- years[-1]
test2 <- data.frame(cbind(years,classes,fbdata))
bubbleplot()
title(main='Winter Pot Survey')

# Create Trawl length vs predicted length freq figure 
bdata <- dat$ont+dat$oot - rept$ent
fbdata <- as.vector(bdata)
classes <- rep(1:6,dat$nyt)
years <-1
for (i in 1:dat$nyt){
years <- c(years,rep(year[dat$it[i]],6))
}
years <- years[-1]
test2 <- data.frame(cbind(years,classes,fbdata))
bubbleplot()
title(main='Trawl Survey')

# Create Pot survey length vs predicted length freq figure 
bdata <- dat$onp+dat$oop - (rept$enp + rept$eop) 
fbdata <- as.vector(bdata)
classes <- rep(1:6,dat$nyp)
years <-1
for (i in 1:dat$nyp){
years <- c(years,rep(year[dat$ip[i]],6))
}
years <- years[-1]
test2 <- data.frame(cbind(years,classes,fbdata))
bubbleplot()
title(main='Summer Pot Survey')

# Create Observer survey length vs predicted length freq figure 
bdata <- dat$ono+dat$ooo - (rept$eno + rept$eoo) 
fbdata <- as.vector(bdata)
classes <- rep(1:6,dat$nyo)
years <-1
for (i in 1:dat$nyo){
years <- c(years,rep(year[dat$io[i]],6))
}
years <- years[-1]
test2 <- data.frame(cbind(years,classes,fbdata))
bubbleplot()
title(main='Observer Survey')

mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)



##############################################################################################################
#  Figure for length stack graph
##############################################################################################################

par(mfrow=c(4,1), mar=c(2,2,2,1),oma=c(3,1,2,1))
color <- c('white', 'gray', 'dimgray', 'darkgray','black','gray90') 
# Create commercial length freq figure 
bdata <- dat$onc+dat$ooc 
barplot(bdata[,-dat$ny],names.arg = year[-dat$ny], col =color)
title(main='Commercial Harvest')

mdata <- matrix(0,6,dat$ny-1)
# Create Winter length freq figure 
bdata <- dat$onw+dat$oow
for (j in 1:dat$nyw) {
	mdata[,dat$iw[j]] <- bdata[,j]
	}
barplot(mdata,names.arg = year[-dat$ny], col =color)	
title(main='Winter Pot Survey')

mdata <- matrix(0,6,dat$ny-1)
# Create Summer Trawl length freq figure 
bdata <- bdata <- dat$ont+dat$oot 
for (j in 1:dat$nyt) {
	mdata[,dat$it[j]] <- bdata[,j]
	}
barplot(mdata,names.arg = year[-dat$ny], col =color)	
title(main='Summter Trawl Survey')

mdata <- matrix(0,6,dat$ny-1)
# Create Summer Pot length freq figure 
bdata <- bdata <- dat$onp+dat$oop 
for (j in 1:dat$nyp) {
	mdata[,dat$ip[j]] <- bdata[,j]
	}
barplot(mdata,names.arg = year[-dat$ny], col=color, legend = c(1,2,3,4,5,6))	
title(main='Summter Pot Survey')

mtext("1: 74-83, 2: 84-93, 3: 94-103, 4: 104-113, 5: 114-123, 6: >124", side=1,line=1.5, outer=T)



















