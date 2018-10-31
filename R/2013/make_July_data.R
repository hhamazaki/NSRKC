#############################################################################################
#   make_July_data.r
#   This program builds a dataset for July schedule model
#   The difference is July model starts July 01 1976 as first model year (year 1), and 
#   winter harvest of Feb 01 1977 is considered occuring in model year 1. 
#   Thus, we shift winter fishery year by 1 year, and add 0 to the final model year's harvest
#############################################################################################

Julydata <- function(dat)
{
# Winter commercial starts 1976 winter
dat$twc <- dat$twc[-1]
dat$tws <- dat$tws[-1]
dat$twst <- dat$twst[-1]
# Add zero to the last year of fishery
dat$twc[dat$ny] <- 0
dat$tws[dat$ny] <- 0
dat$twst[dat$ny] <- 0

# Winter Pot survey is also shifted 1 year
dat$iw <- dat$iw-1
return(dat)
}










