########################################################################
## Interpolation of of event data (five year periods) 
## fifth/single years by Beers Modified six-term formula
## (Siegel and Swanson, 2004, p. 729)
## This formula applies some smoothing to the interpolant, which is
## recommended for time series of events with some dynamic
## R implementation by Thomas Buettner (21 Oct. 2015)
########################################################################

# save working directory
wd <- getwd()

#set working directory:
setwd("v:/R/Functions/interpolation/BeersModifiedSplit")
                                                      
# input
fn <- "Births"
ifn <- paste(fn, "5.csv", sep = "")     # file name input

# output
ofn <- paste(fn, "1.csv", sep = "")     # file name output

tp <- read.csv(ifn, header = TRUE, row.names = 1)

## interpolation period
YEAR1 <- as.numeric(substr(rownames(tp),1,4))
YEAR2 <- as.numeric(substr(rownames(tp),6,10))
YEAR <- 0.5 + (YEAR1 + YEAR2)/2

FYEAR <- min(YEAR1) #1950
LYEAR <- max(YEAR2) #2100
LYEAR1 <- LYEAR - 1

## dimensions
NCOL <- dim(tp)[2]   ## countries or trajectories
NAG5 <- dim(tp)[1]   ## age or time

# number of single year or single age groups
NAG1 <- NAG5 * 5

## Beers Modified Split
bc <- c(
  0.3332, -0.1938,  0.0702, -0.0118,  0.0022 ,
  0.2569, -0.0753,  0.0205, -0.0027,  0.0006 ,
  0.1903,  0.0216, -0.0146,  0.0032, -0.0005 ,
  0.1334,  0.0969, -0.0351,  0.0059, -0.0011 ,
  0.0862,  0.1506, -0.0410,  0.0054, -0.0012 ,
  
  0.0486,  0.1831, -0.0329,  0.0021, -0.0009 ,
  0.0203,  0.1955, -0.0123, -0.0031, -0.0004 ,
  0.0008,  0.1893,  0.0193, -0.0097,  0.0003 ,
 -0.0108,  0.1677,  0.0577, -0.0153,  0.0007 ,
 -0.0159,  0.1354,  0.0972, -0.0170,  0.0003 ,
  
 -0.0160,  0.0973,  0.1321, -0.0121, -0.0013 ,
 -0.0129,  0.0590,  0.1564,  0.0018, -0.0043 ,
 -0.0085,  0.0260,  0.1650,  0.0260, -0.0085 ,
 -0.0043,  0.0018,  0.1564,  0.0590, -0.0129 ,
 -0.0013, -0.0121,  0.1321,  0.0973, -0.0160 ,
  
  0.0003, -0.0170,  0.0972,  0.1354, -0.0159 ,
  0.0007, -0.0153,  0.0577,  0.1677, -0.0108 ,
  0.0003, -0.0097,  0.0193,  0.1893,  0.0008 ,
 -0.0004, -0.0031, -0.0123,  0.1955,  0.0203 ,
 -0.0009,  0.0021, -0.0329,  0.1831,  0.0486 ,
  
 -0.0012,  0.0054, -0.0410,  0.1506,  0.0862 ,
 -0.0011,  0.0059, -0.0351,  0.0969,  0.1334 ,
 -0.0005,  0.0032, -0.0146,  0.0216,  0.1903 ,
  0.0006, -0.0027,  0.0205, -0.0753,  0.2569 ,
  0.0022, -0.0118,  0.0702, -0.1938,  0.3332 
)
## standard format for Beers coefficients
bm <- matrix(bc,25,5,byrow = T)

## Age vector
a <- seq(0, (NAG5 - 1)*5, by = 5)

# number of middle panels
MP <- NAG5 - 5 + 1

## creating a beers coefficient matrix for 18 5-year age groups
bcm <- array(0, dim = c(NAG1,NAG5))

## insert first two panels
bcm[1:10,1:5] <- bm[1:10,]

for (i in (1:MP))
{
  # calculate the slices and add middle panels accordingly
  bcm[((i + 1)*5 + 1):((i + 2)*5), i:(i + 4)] <- bm[11:15,]
}

## insert last two panels
bcm[((NAG5 - 2)*5 + 1):(NAG5*5),MP:(MP + 4)] <- bm[16:25,]

options(max.print = 10000000) 

pop <- bcm %*% as.matrix(tp)

## write un-shifted results
## write.csv(pops, ofn)

## shifting the interpolant by half a year forward
## In general, each interval is halved and the the halves are recombined to form a calendar year.
## The first missing half year is assumed to be equal the half of the first interval.
## This means simply that the first year is equal to the first (shifted) interval
## The implementation used vector arithmetic by constructing a parallel data array 
## with the first element equal to the original first interval, 
## and filling the rest with the original data up to n-1 elements
## Adding the two data structures and dividing the result by two 
## produces the shifted time reference.

pop2 <- array(0, dim = c(NAG1,NCOL))
pop2[1,] <- pop[1,]
pop2[2:NAG1,] <- pop[1:(NAG1 - 1),]
pops <- (pop + pop2)*0.5
## write shifted results
write.csv(pops, ofn)