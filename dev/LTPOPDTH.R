
# Author: tim
###############################################################################

# pars
sex    <- "m"
IMR    <- 0.1 # or NA if unknown

# separation factors (what is this?)

code   <- 2
# for empirical separation factors?
Age0   <- NA   # looks like these are ax values for ages 0 and 1-4...
Age1_4 <- NA

SRB    <- 1.05 # used for both-sex lifetables only

# primary inputs are pop and exposures in abridged ages. Abridged lifetables
# used. Not a fan: they use 2.5 ax... leave arg open for axmethod, or to 
# abridge a single age table the way the HMD does?
# can use calcAgeAbr() plus groupAges() to move single to abridged.
# 


Exposures <- c(100958,466275,624134,559559,446736,370653,301862,249409,
247473,223014,172260,149338,127242,105715,79614,53660,
31021,16805,8000,4000,2000,1000)

Deaths <- c(8674,1592,618,411,755,1098,1100,1357,
1335,3257,2200,4023,2167,4578,2956,4212,
2887,2351,1500,900,500,300)

# lower age bounds
Age    <- c(0, 1, seq(5, 100, by = 5))
AgeInt <- c(diff(Age), NA)

# some calcs
nMx <- Deaths / Exposures

# 
source("/home/tim/git/DemoTools/dev/nAx.R")
nAx <- nMx2nAxCDabr(nMx, AgeInt = AgeInt, IMR = NA, sex = "m", region = "n")

mxax2qx <- function(nax, nMx, AgeInt, closeout = TRUE, IMR){
	
	qx <- (AgeInt * nMx) / (1 + (AgeInt - nax) * nMx)
	if (closeout){
		qx[length(qx)] <- 1
	}
	if (!missing(IMR) & !is.na(IMR)){
		qx[1] <- IMR
	}
	qx
}
nqx <- mxax2qx(nAx, nMx, AgeInt, TRUE, .1)

qx2lx <- function(nqx,radix=1e5){
	
	radix * cumprod(c(1, 1 - nqx[-length(nqx)]))
}

lx <- qx2lx(nqx)

lx2dx <- function(lx){
	diff(-c(lx,0))
}
ndx <- lx2dx(lx)
plot(ndx)
