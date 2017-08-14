
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


qx2mx <- function(nqx, ax, AgeInt = inferAgeIntAbr(vec = nqx)) {
	nqx / (AgeInt - (AgeInt - ax) * nqx)      }

mx2qx <- function(mx, ax, AgeInt = inferAgeIntAbr(vec = mx)) {
	(AgeInt * mx) / (1.0 + ((AgeInt - ax) * mx)) 
}

qxmx2ax <- function(nqx, nmx, AgeInt = inferAgeIntAbr(vec = nqx)){
	1 / nmx - AgeInt / nqx + AgeInt
}

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




# some calcs


LTabr <- function(
		Deaths, 
		Exposures, 
		nMx,
		Age,
		AgeInt = inferAgeIntAbr(Age = Age), 
		axmethod = "pas", 
		Sex = "m", 
		region = "w",
		SRB = 1.05, 
		IMR){
	
	if (missing(Deaths)){
		Deaths    <- NULL
	}
	if (missing(Exposures)){
		Exposures <- NULL
	}
	
	axmethod      <- tolower(axmethod)
	Sex           <- tolower(Sex)
	region        <- tolower(region)
	
	if (missing(nMx)){
		nMx       <- Deaths / Exposures
    }
	if (missing(IMR)){
		IMR       <- nMx[1]
	}
	
	# take care of ax first, two ways presently
	if (axmethod == "pas"){
	    nAx       <- axPAS(
				      nMx, 
				      AgeInt = AgeInt, 
				      IMR = IMR, 
				      Sex = Sex, 
				      region = region)
    }
	if (axmethod == "un"){
		# UN method just CD west for now, so no region arg
		nAx       <- axUN(nMx = nMx, 
				          AgeInt = AgeInt, 
						  Sex = Sex)
	}
    # these are the same always. Only thing different is ax method.
	nqx <- mxax2qx(nAx, nMx = nMx, AgeInt = AgeInt, closeout = TRUE, IMR = IMR)
	lx  <- qx2lx(nqx)
	ndx <- lx2dx(lx)
	nLx <- lxdxax2Lx(lx, ndx, nAx, AgeInt)
	Tx  <- Lx2Tx(nLx)
	ex  <- Tx / lx
	
	# output is an unrounded, unsmoothed lifetable
	out <- data.frame(Deaths = Deaths, 
			   Exposures = Exposures,
			   Age = Age,
			   AgeInt = AgeInt,
			   nMx = nMx,
			   nAx = nAx,
			   nqx = nqx,
			   lx = lx,
			   ndx = ndx,
			   nLx = nLx,
			   Tx = Tx,
			   ex = ex)
	return(out)
}




# trial code from PAS LTPOPDTH, North, Males, IMR = .1
Exposures <- c(100958,466275,624134,559559,446736,370653,301862,249409,
		247473,223014,172260,149338,127242,105715,79614,53660,
		31021,16805,8000,4000,2000,1000)

Deaths <- c(8674,1592,618,411,755,1098,1100,1357,
		1335,3257,2200,4023,2167,4578,2956,4212,
		2887,2351,1500,900,500,300)
# lower age bounds
Age    <- c(0, 1, seq(5, 100, by = 5))
AgeInt <- c(diff(Age), NA)



PASLT <- LTabr(Deaths = Deaths, 
		Exposures = Exposures, 
		Age = Age,
		AgeInt =AgeInt,
		axmethod = "PAS",
		IMR = .1,
		region = "N",
		Sex = "M")
# de facto unit test. The unsmoothed output from PAS
# spreadsheet (allow for rounding error in last decimal)
excheck <- c(56.31,61.53,58.35,53.63, 
		48.81,44.21,39.83,35.52,
		31.43,27.22,24.09,20.52,
		18.12,14.51,12.42,9.45,7.85, 
		6.09,4.95,4.28,3.85,3.33)
stopifnot(max(abs(round(PASLT$ex,2) - excheck)) <= .01)

# example data from United Nations Model Life Tables -
# Males Latin American Pattern

Mx <- c(.23669,.04672,.00982,.00511,.00697,.01036,.01169,
.01332,.01528,.01757,.02092,.02517,.03225,.04241,.06056,
.08574,.11840,.16226,.23745)
qx <- c(.20429,.16631,.04790,.02522,.03427,.05051,.05679,.06449,
  .07363,.08418,.09948,.11849,.14939,.19205,.26327,.35208,
  .45210,.56382,1
 )
ax <- c(0.330,1.352,2.500,2.500,2.633,2.586,2.528,2.528,
 2.526,2.529,2.531,2.538,2.542,2.543,2.520,2.461,2.386,2.295,4.211)

