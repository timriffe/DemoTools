
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



# 
source("/home/tim/git/DemoTools/dev/nAx.R")


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


qx2lx <- function(nqx,radix=1e5){
	
	radix * cumprod(c(1, 1 - nqx[-length(nqx)]))
}

lx2dx <- function(lx){
	diff(-c(lx,0))
}
# dx <- ndx; ax <- nAx
lxdxax2Lx <- function(lx,dx,ax,AgeInt){
	N                   <- length(lx)
	Lx                  <- rep(0,N)
	Lx[1:(N - 1)]       <- AgeInt[1:(N - 1)] * lx[2:N] + ax[1:(N - 1)] * dx[1:(N - 1)]
	Lx[N]		        <- lx[N] * ax[N]
	Lx
}
Lx2Tx <- function(Lx){
	rev(cumsum(rev(nLx)))
}
# trial code:

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
nMx    <- Deaths / Exposures

nAx <- axPAS(nMx, AgeInt = AgeInt, IMR = NA, sex = "m", region = "n")

# these are the same always. Only thing different is ax method.
nqx <- mxax2qx(nAx, nMx, AgeInt, closeout = TRUE, IMR = .1)
lx  <- qx2lx(nqx)
ndx <- lx2dx(lx)
nLx <- lxdxax2Lx(lx, ndx, nAx, AgeInt)
Tx  <- Lx2Tx(nLx)
ex  <- Tx / lx

LTabr <- function(
		Deaths, 
		Exposures, 
		Age,
		AgeInt = inferAgeIntAbr(Age = Age), 
		axmethod = "pas", 
		sex = "m", 
		region = "w",
		SRB = 1.05, 
		IMR = Deaths[1] / Exposures[1], 
		...){
	
	axmethod <- tolower(axmethod)
	sex      <- tolower(sex)
	region   <- tolower(region)
	
	nMx      <- Deaths / Exposures
	
	if (axmethod = "pas"){
	    nAx <- axPAS(nMx, AgeInt = AgeInt, IMR = IMR, sex = sex, region = region)
    }
	if (axmethod = "un"){
		
	}
# these are the same always. Only thing different is ax method.
	nqx <- mxax2qx(nAx, nMx = nMx, AgeInt = AgeInt, closeout = TRUE, IMR = IMR)
	lx  <- qx2lx(nqx)
	ndx <- lx2dx(lx)
	nLx <- lxdxax2Lx(lx, ndx, nAx, AgeInt)
	Tx  <- Lx2Tx(nLx)
	ex  <- Tx / lx
	
	data.frame(Deaths = Deaths, 
			   Exposures = Exposures,
			   Age)
}

plot(ndx)
