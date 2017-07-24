########################################################################
## Interpolation of five year age groups into single year age groups by 
## Sprague's fifth-difference  osculatory formula
## (Siegel and Swanson, 2004, p. 727)
## R implementation by Thomas Buettner (21 Oct. 2015)
## expanded by Patrick Gerland (18 feb. 2016): 
## - added monotonic decline for age 90+ (and to avoid negative values at upper ages)
## - added cohort interpolation between pivotal years ending in 0 and 5
## - implemented exponential growth rate cohort interpolation instead of linear growth (WPP default)
## - enforced consistency between annually interpolated total using Beers 
## and interpolated pop1x1 for intermediate years between pivotal years (use only age distribution)
########################################################################

## used for pchip interpolation
if(!require(signal)){
    install.packages("signal")
    library(signal)
}

# save working directory
wd <- getwd()

# set working directory:
setwd("C:/Users/Patrick/Dropbox/R/interpolation/wpp/")
ifn <- "/home/tim/git/DemoTools/DemoTools/dev/Data/IND_Pop2015B1_1.7.csv"
# read data to be interpolated
# input data are expected to have age as rows, and countries as columns. 
# The first column has the row labels (=age), the first row holds column labels (countrIDs, year...) 
# it is assumed that input data have a last, open-ended age group which must not be interpolated.
#
# file name convention is shown below. 
# In the example, the original inpout file is named WPP2012Pop5.csv, or *5.csv
# The script produces two outputs:
# *1.csv -> Interpolated data tabulated like to the input data (over time and by age)
# *1c.csv -> Interpolated data tabulated like to the input data (over time and by age, and then by cohort)

# input

fn    <- "IND_Pop2015B" # annually interpolated using Beers in column
ifn   <- paste(fn, "1_1.7.csv", sep = "")

# output
ofn   <- paste(fn, "1x1_1.7.csv", sep = "")
ofncl <- paste(fn, "1x1-cohort-linear-interpol.csv", sep = "")  ## cohort interpolated
ofnc  <- paste(fn, "1x1-cohort-exponential-interpol.csv", sep = "")  ## cohort interpolated

## read annually interpolated data (transposed by age in rows and time in column)
tp    <- t(read.csv(ifn, header = TRUE, row.names = 1, check.names = FALSE))

Age5 <- as.numeric(row.names(tp))
FAGE <- Age5[1]
LAGE <- Age5[length(Age5)]
Age1 <- seq(FAGE, LAGE, by = 1)
# Age1

## interpolation period

YEAR   <- as.numeric(substr(colnames(tp),1,5))
FYEAR  <- min(YEAR) #1950
LYEAR  <- max(YEAR) #2100
LYEAR1 <- LYEAR

## dimensions
NCOL   <- dim(tp)[2]
NAG5   <- dim(tp)[1] 

# number of single year age groups (closed age groups only)
NAG1   <- NAG5 * 5 - 5 

## Sprague's Split
bc <- c(
  0.3616, -0.2768,  0.1488, -0.0336,  0.0000, 
  0.2640, -0.0960,  0.0400, -0.0080,  0.0000, 
  0.1840,  0.0400, -0.0320,  0.0080,  0.0000, 
  0.1200,  0.1360, -0.0720,  0.0160,  0.0000, 
  0.0704,  0.1968, -0.0848,  0.0176,  0.0000, 
  0.0336,  0.2272, -0.0752,  0.0144,  0.0000, 
  0.0080,  0.2320, -0.0480,  0.0080,  0.0000, 
  -0.0080,  0.2160, -0.0080,  0.0000,  0.0000, 
  -0.0160,  0.1840,  0.0400, -0.0080,  0.0000, 
  -0.0176,  0.1408,  0.0912, -0.0144,  0.0000, 
  -0.0128,  0.0848,  0.1504, -0.0240,  0.0016, 
  -0.0016,  0.0144,  0.2224, -0.0416,  0.0064, 
  0.0064, -0.0336,  0.2544, -0.0336,  0.0064, 
  0.0064, -0.0416,  0.2224,  0.0144, -0.0016, 
  0.0016, -0.0240,  0.1504,  0.0848, -0.0128, 
  0.0000, -0.0144,  0.0912,  0.1408, -0.0176, 
  0.0000, -0.0080,  0.0400,  0.1840, -0.0160, 
  0.0000,  0.0000, -0.0080,  0.2160, -0.0080, 
  0.0000,  0.0080, -0.0480,  0.2320,  0.0080, 
  0.0000,  0.0144, -0.0752,  0.2272,  0.0336, 
  0.0000,  0.0176, -0.0848,  0.1968,  0.0704, 
  0.0000,  0.0160, -0.0720,  0.1360,  0.1200, 
  0.0000,  0.0080, -0.0320,  0.0400,  0.1840, 
  0.0000, -0.0080,  0.0400, -0.0960,  0.2640, 
  0.0000, -0.0336,  0.1488, -0.2768,  0.3616
)
## standard format for Sprague coefficients
sm  <- matrix(bc, 25, 5, byrow = TRUE)

# number of middle panels
MP  <- NAG5 - 5 

## creating a Spragues coefficient matrix for 5-year age groups
scm <- array(0, dim = c(NAG1 + 1, NAG5))

## insert first two panels
scm[1:10,1:5] <- sm[1:10,]

for (i in (1:MP)) {
  # calculate the slices and add middle panels accordingly
  scm[((i + 1)*5 + 1):((i + 2)*5), i:(i + 4)] <- sm[11:15,]
}



## insert last two panels
fr <- (NAG5 - 3) * 5 + 1
lr <- (NAG5 - 1) * 5

fc <- MP 
lc <- MP + 4 

#scm[((NAG5 - 3)*5 + 1):((NAG5-1)*5),MP - 1:(MP-1 + 4)] <- sm[16:25,]
scm[fr:lr,fc:lc] <- sm[16:25,]
# last open ended age group
scm[NAG1 + 1,NAG5] <- 1


options(max.print = 10000000) 

pop           <- scm %*% as.matrix(tp)
rownames(pop) <- Age1


## write intermediate results for checking
write.csv(pop, paste(fn, "1x1-step1-sprague.csv", sep = ""))


##############################################################################
# add check and treatment for negative values (usually at highest ages) here
# P. Gerland (18 Feb 2016) -- Note: don't know what/how WPP .VB code is doing for this step 
# but conceptually this approach is robust and effective in imposing a monotonically 
# declining population at upper ages and preventing negative values.
##############################################################################

pop.pchip <- NULL
for (j in 1:NCOL){
	## interpolate on cumulated distribution
	yinterpol <- interp1(c(0, Age5+5), c(0, cumsum(tp[,j])), c(0, Age1+1), 'pchip', extrap = TRUE)
	## decumulate
	yinterpol <- diff(yinterpol)
	pop.pchip <- cbind(pop.pchip, yinterpol)
}
colnames(pop.pchip) <- colnames(pop)

## tail(pop, 12)
## tail(pop.pchip, 12)

## combine Sprague interpolation up to age 90 with pchip interpolation for age 90+

## 90-94
# tp[19,] # Pop 90-94
# pop[(91:95),] # Sprague
# pop.pchip[(91:95),] # pchip

## Pivotal age 90
## impute pchip interpolation as default
pop.combined     <- pop.pchip[91,]
## substitute Sprague interpolation if > pchip for better belnding of the two series
pop.combined[pop[91,] > pop.pchip[91,]] <- pop[91,pop[91,] > pop.pchip[91,]]
## combine age 91-94 from pchip
pop.combined     <- rbind(pop.combined, pop.pchip[(92:95),])
## sum of rows
pop.combined.sum <- colSums(pop.combined, na.rm = TRUE)
## adjust back on initial pop 5x5 for age 90-94
## proportional distribution
pop.combined[is.na(pop.combined)==TRUE] <- 0
prop             <- prop.table(as.matrix(pop.combined), margin=2)
rownames(prop)   <- c(91:95)
df               <- as.data.frame(prop)
## vector of initial pop 5x5 for age 90-94
v                <- as.vector(tp[19,])
pop.combined     <- data.frame(mapply(`*`,df,v))

## append the remaining of the age groups (except last open age)
## 95-99 onward
# tp[(20:(NAG5-1)),] ## pop5
# pop[(96:NAG1),]
# pop.pchip[(96:NAG1),]
colnames(pop.combined) <- colnames(pop.pchip)
pop.combined <- rbind(pop.combined, pop.pchip[(96:NAG1),])
pop.combined <- rbind(pop.combined, pop[(NAG1+1),])

## append Sprague interpolation before age 90
pop.combined <- rbind(pop[(1:90),], pop.combined)

## deal with negative values if applicable (but in principle should not be happening)
pop.combined[pop.combined<0] <- 0
rownames(pop.combined) <- Age1

## write intermediate results for checking
write.csv(pop.combined, paste(fn, "1x1-step2-pchipGE90.csv", sep = ""))



## extra step: P. Gerland - cohort interpolation on WPP pivotal years ending in 0 and 5 for age 0-94 in year t to age 5-99 in year t+5
## use exponential growth rate (r=ln(pt2/pt1)/5 with ptn = pt1*(exp(r*n)) ) instead of linear interpolation (WPP default which causes convex hump by age over time at age 80+)
## keep pop1x1 interpolation for intermediate years for cohorts < age 5 in year t+5 not available in year t
## keep pop1x1 interpolation for intermediate years for cohorts age 95+ in year t not available in year t+5

## get pivotal years
index    <- seq(from = 1, to = NCOL, by = 5)
n.index  <- length(index)
data     <- pop.combined[,index]

## get Sprague interpolated pop1x1 for intermediate years (used to impute pop for cohorts that cannot be interpolated)
sprague1 <- pop.combined[,(index+1)[1:(n.index-1)]]
sprague2 <- pop.combined[,(index+2)[1:(n.index-1)]]
sprague3 <- pop.combined[,(index+3)[1:(n.index-1)]]
sprague4 <- pop.combined[,(index+4)[1:(n.index-1)]]

## pop1 in year t
data1    <- data[1:(NAG1-5),]
## head(data1, 10)
## tail(data1, 10)

## pop2 in year t+5
data2    <- data[6:NAG1,]
## head(data2, 10)
## tail(data2, 10)

## growth rate between yeart (t, t+5)
r.linear <- (data2[,2:length(index)] / data1[,1:(length(index)-1)]) / 5
r        <- log(data2[,2:length(index)] / data1[,1:(length(index)-1)]) / 5
## head(r)

## cohort interpolated pop in year t+1 to t+4
popt1.linear <- data1[,1:(n.index-1)] + 1 * (data2[,2:length(index)] - data1[,1:(length(index)-1)]) / 5
popt2.linear <- data1[,1:(n.index-1)] + 2 * (data2[,2:length(index)] - data1[,1:(length(index)-1)]) / 5
popt3.linear <- data1[,1:(n.index-1)] + 3 * (data2[,2:length(index)] - data1[,1:(length(index)-1)]) / 5
popt4.linear <- data1[,1:(n.index-1)] + 4 * (data2[,2:length(index)] - data1[,1:(length(index)-1)]) / 5

popt1        <- data1[,1:(n.index-1)] * exp(r*1)
popt2        <- data1[,1:(n.index-1)] * exp(r*2)
popt3        <- data1[,1:(n.index-1)] * exp(r*3)
popt4        <- data1[,1:(n.index-1)] * exp(r*4)
## head(popt1)

## combined dataset
## linear
pop1x1.linear                                       <- pop.combined  ## initial pop1x1 from Sprague
## t+1 substitute with popt1
pop1x1.linear[2:(NAG1-4), (index+1)[1:(n.index-1)]] <- popt1.linear
## t+2 substitute with popt2
pop1x1.linear[3:(NAG1-3), (index+2)[1:(n.index-1)]] <- popt2.linear
## t+3 substitute with popt3
pop1x1.linear[4:(NAG1-2), (index+3)[1:(n.index-1)]] <- popt3.linear
## t+4 substitute with popt4
pop1x1.linear[5:(NAG1-1), (index+4)[1:(n.index-1)]] <- popt4.linear

## exponential
pop1x1 <- pop.combined  ## initial pop1x1 from Sprague
## t+1 substitute with popt1
pop1x1[2:(NAG1-4), (index+1)[1:(n.index-1)]] <- popt1 
## t+2 substitute with popt2
pop1x1[3:(NAG1-3), (index+2)[1:(n.index-1)]] <- popt2
## t+3 substitute with popt3
pop1x1[4:(NAG1-2), (index+3)[1:(n.index-1)]] <- popt3
## t+4 substitute with popt4
pop1x1[5:(NAG1-1), (index+4)[1:(n.index-1)]] <- popt4


## check
## round(head(pop1x1, 101)[,(NCOL-6):NCOL],3)


## insure consistency on annual total (interpolated using Beers) by sex for intermediate years

## linear cohort interpolation
## proportional distribution
prop                      <- prop.table(as.matrix(pop1x1.linear), margin=2)
rownames(prop)            <- rownames(pop1x1.linear)
df                        <- as.data.frame(prop)
## vector of total pop by sex to use to enforce overall consistency by year (annually interpolated Beers total by sex)
total                     <- apply(tp, 2, sum, na.rm=TRUE)

pop1x1.adjusted           <- data.frame(mapply(`*`,df,total))
## combine back with IDs
rownames(pop1x1.adjusted) <- rownames(pop1x1)
colnames(pop1x1.adjusted) <- colnames(pop1x1)
## Round to significant digit
output.wide               <- round(pop1x1.adjusted, 3)

## write results
write.csv(pop1x1.adjusted, ofncl)
## write results (Rounded to significant digit)
write.csv(output.wide, ofncl)


## exponential cohort interpolation
## proportional distribution
prop                      <- prop.table(as.matrix(pop1x1), margin=2)
rownames(prop)            <- rownames(pop1x1)
df                        <- as.data.frame(prop)
## vector of total pop by sex to use to enforce overall consistency by year (annually interpolated Beers total by sex)
total                     <- apply(tp, 2, sum, na.rm=TRUE)

pop1x1.adjusted           <- data.frame(mapply(`*`,df,total))
## combine back with IDs
rownames(pop1x1.adjusted) <- rownames(pop1x1)
colnames(pop1x1.adjusted) <- colnames(pop1x1)
## Round to significant digit
output.wide <- round(pop1x1.adjusted, 3)

## write results
write.csv(pop1x1.adjusted, ofnc)
## write results (Rounded to significant digit)
write.csv(output.wide, ofnc)


colSums(tp) - 
colSums(pop)



tpmini      <- tp[, 1, drop = FALSE]
tpmini[6, ] <- colSums(tpmini[6:nrow(tp), , drop = FALSE])
tpmini      <- tpmini[1:6, , drop = FALSE]
colSums(spragueSimple(tpmini)) - colSums(tpmini)


colSums(spragueSimple(tpmini)) - colSums(tpmini)
colSums(spragueSimple(tp)) - colSums(tp)
popmat <- tp

image(scm[,1:5]- scm2[,1:5])

#head(scm,10) - 
#head(scm2,10)
#scm[11, ]
#scm2[10, ]