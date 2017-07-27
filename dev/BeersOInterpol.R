########################################################################
## Interpolation of quinquennial stock data into fifth/single years by 
## Beers Ordinary fifth-difference  osculatory formula
## (Siegel and Swanson, 2004, p. 728)
## R implementation by Thomas Buettner (21 Oct. 2015)
########################################################################
## PG: used in WPP to interpolate total pop by sex over time, and pop for each age group over time in crosssectional way (alternative to ln(p2/p1))/(t2-t1)) -- smoother trend -- preserves input value in pivotal year ending in 0 and 5

# save working directory
wd <- getwd()

#set working directory:
setwd("C:/Users/Patrick/Dropbox/R/interpolation/wpp/")

# input 
## rows used for 5-year spaced stock data and columns for age groups or countries/areas or trajectories
fn <- "IND_Pop2015B" # pop by 5-year age groups in column

## input data by 5-year time intervals
ifn <- paste(fn, "5_1.7.csv", sep = "")

# output, annual mid-year
ofn <- paste(fn, "1_1.7.csv", sep = "")

#output, annual beginning of year
ofn1 <- paste(fn, "1_1.1.csv", sep = "")

tp <- read.csv(ifn, header = TRUE, row.names = 1, check.names = FALSE)

YEAR <- as.numeric(substr(rownames(tp),1,4))
FYEAR <- min(YEAR) #1950
LYEAR <- max(YEAR) #2100


## dimensions
NY <- 5                        # Number of years in period
NCOL <- dim(tp)[2]             # number of data columns/number of countries/time series
NYEAR5 <- dim(tp)[1]           # number of five year period
NYEAR1 <- (NYEAR5 - 1) * 5 + 1 # number of single years

## Beers Ordinary Interpol
bc <- c(
 1.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
 0.6667,  0.4969, -0.1426, -0.1006,  0.1079, -0.0283,
 0.4072,  0.8344, -0.2336, -0.0976,  0.1224, -0.0328,
 0.2148,  1.0204, -0.2456, -0.0536,  0.0884, -0.0244,
 0.0819,  1.0689, -0.1666, -0.0126,  0.0399, -0.0115,
 0.0000,  1.0000,  0.0000,  0.0000,  0.0000,  0.0000,
-0.0404,  0.8404,  0.2344, -0.0216, -0.0196,  0.0068,
-0.0497,  0.6229,  0.5014, -0.0646, -0.0181,  0.0081,
-0.0389,  0.3849,  0.7534, -0.1006, -0.0041,  0.0053,
-0.0191,  0.1659,  0.9354, -0.0906,  0.0069,  0.0015,
 0.0000,  0.0000,  1.0000,  0.0000,  0.0000,  0.0000,
 0.0117, -0.0921,  0.9234,  0.1854, -0.0311,  0.0027,
 0.0137, -0.1101,  0.7194,  0.4454, -0.0771,  0.0087,
 0.0087, -0.0771,  0.4454,  0.7194, -0.1101,  0.0137,
 0.0027, -0.0311,  0.1854,  0.9234, -0.0921,  0.0117,
 0.0000,  0.0000,  0.0000 , 1.0000,  0.0000,  0.0000,
 0.0015,  0.0069, -0.0906 , 0.9354,  0.1659, -0.0191,
 0.0053, -0.0041, -0.1006 , 0.7534,  0.3849, -0.0389,
 0.0081, -0.0181, -0.0646 , 0.5014,  0.6229, -0.0497,
 0.0068, -0.0196, -0.0216 , 0.2344,  0.8404, -0.0404,
 0.0000,  0.0000,  0.0000,  0.0000,  1.0000,  0.0000,
-0.0115,  0.0399, -0.0126, -0.1666,  1.0689,  0.0819,
-0.0244,  0.0884, -0.0536, -0.2456,  1.0204,  0.2148,
-0.0328,  0.1224, -0.0976, -0.2336,  0.8344,  0.4072,
-0.0283,  0.1079, -0.1006, -0.1426,  0.4969,  0.6667,
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  1.0000
)
## Format for Beers ordinary coefficients 26 x 6 fixed
bm <- matrix(bc,26,6,byrow = T)

# number of middle panels
MP <- NYEAR5 - 5

## creating a beers coefficient matrix for single years
bcm <- array(0, dim = c(NYEAR1,NYEAR5))

## inserting first two panels
bcm[1:10,1:6] <- bm[1:10,]

# inserting middle panels
for (i in (1:MP))
{
  # calculate the slices and add middle panels accordingly
  bcm[((i + 1)*NY + 1):((i + 2)*NY), i:(i + NY)] <- bm[11:15,]
}

## insert last two panels
bcm[((NYEAR5 - 3)*NY + 1):(NYEAR1),MP:(MP + NY)]<- bm[16:26,]

options(max.print = 10000000) 

pop <- bcm %*% as.matrix(tp)

colnames(pop) <- colnames(tp)
rownames(pop) <- c(FYEAR:LYEAR)

write.csv(pop, ofn, row.names = TRUE)



############################
# # Moving the interpolated figures back by half a year, from mid-year to beginning of year
# # from 1.7(t) back to 1.1(t)
pt1 <- array(0, dim = dim(pop))
pt2 <- array(0, dim = dim(pop))
pop11 <- array(0, dim = dim(pop))

pt1 <- pop[1:(nrow(pop)-1),]
pt2 <- pop[2:nrow(pop),]
r2 <- log(pt2/pt1)*0.5
pop11[1:(nrow(pop)-1),] <- pt1 * exp(-r2)
pop11[nrow(pop),] <- pop[nrow(pop),] * exp(+r2[nrow(r2)])

colnames(pop11) <- colnames(tp)
rownames(pop11) <- c(FYEAR:LYEAR)

###############################
write.csv(pop11, ofn1, row.names = TRUE)
