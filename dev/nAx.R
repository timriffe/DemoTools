
# Author: tim
###############################################################################

# TR: lots of approaches to a(0). Sometimes called 'separation factors'. These ought to be 
# collected, organized, and standardized here.

# these are the PAS versions:
# this is called a separation factor in the spreadsheet?
# separate estimate of IMR optional
geta0CD <- function(M0, IMR, sex = "M", region = "W"){
	# sex can be "m", "f", or "b"
	# region can be "n","e","s","w",or
	sex    <- tolower(sex)
	region <- tolower(region)
	
	Age0 <- matrix(c(0.33, 	0.35, 	0.3500, 
					0.33, 	0.35, 	0.3500, 
					0.29, 	0.31, 	0.3100,
					0.33, 	0.35, 	0.3500 
			), ncol = 3, byrow = TRUE,
			dimnames = list(c("w", "n", "e", "s"), c("m","f","b")))
	
	# IMR optional here
	if (missing(IMR) | is.na(IMR)){
		IMR <- M0
	}
	ifelse(IMR > .1, Age0[region, sex], .33)
}
# separate estimate of IMR optional
# TR: I think it's funny that a1-4 doesn't depend at all on m1-4
geta1_4CD <- function(M0, IMR, sex = "M", region = "W"){
	# sex can be "m", "f", or "b"
	# region can be "n","e","s","w",or
	sex    <- tolower(sex)
	region <- tolower(region)
	Age1_4 <- matrix(c(1.352, 	1.361, 	1.3610, 
					1.558, 	1.570, 	1.5700, 
					1.313, 	1.324, 	1.3240, 
					1.240, 	1.239, 	1.2390
			), ncol = 3, byrow = TRUE,
			dimnames = list(c("w", "n", "e", "s"), c("m","f","b")))
	if (!missing(IMR) & !is.na(IMR)){
		q0 <- IMR
	} else {
		a0 <- geta0CD(M0, sex = sex, region = region)
		q0 <- (M0 * a0) / 
				(1 + (AgeInt[1] - a0 * M0))
	}
	ifelse(q0 > .1, Age1_4[region, sex], 1.5577)
}

# so need args for a0 and a1_4
# make a single abridged ax function that does this stuff
# this is a terrible name, hate it. Not elegant
nMx2nAxCDabr <- function(nMx, AgeInt, IMR, sex = "M", region = "W", OAG = TRUE){
	# sex can be "m", "f", or "b"
	# region can be "n","e","s","w",or
	sex    <- tolower(sex)
	region <- tolower(region)
	
	N      <- length(nMx)
	a0     <- rep(2.5, N)
	
	if (missing(IMR) | is.na(IMR)){
		IMR <- nMx[1]
	}
	a0[1]  <- geta0CD(M0 = nMx[1], IMR = IMR, sex = sex, region = region)
	a0[2]  <- geta1_4CD(M0 = nMx[1], IMR = IMR, sex = sex, region = region)
	a0[N]  <- ifelse(!is.na(AgeInt[N]) & AgeInt[N] == 5, 2.5, 1 / nMx[N])
	a0
}