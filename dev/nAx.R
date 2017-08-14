
# Author: tim
###############################################################################

# TR: lots of approaches to a(0). Sometimes called 'separation factors'. These ought to be 
# collected, organized, and standardized here.

# these are the PAS versions:
# this is called a separation factor in the spreadsheet?
# separate estimate of IMR optional
geta0CD <- function(M0, IMR, Sex = "M", region = "W"){
	# sex can be "m", "f", or "b"
	# region can be "n","e","s","w",or
	Sex    <- tolower(Sex)
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
	ifelse(IMR > .1, Age0[region, Sex], .33)
}
# separate estimate of IMR optional
# TR: I think it's funny that a1-4 doesn't depend at all on m1-4
geta1_4CD <- function(M0, IMR, Sex = "M", region = "W"){
	# sex can be "m", "f", or "b"
	# region can be "n","e","s","w",or
	Sex    <- tolower(Sex)
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
		a0 <- geta0CD(M0, Sex = Sex, region = region)
		q0 <- (M0 * a0) / 
				(1 + (AgeInt[1] - a0 * M0))
	}
	ifelse(q0 > .1, Age1_4[region, Sex], 1.5577)
}

# so need args for a0 and a1_4
# make a single abridged ax function that does this stuff
# this is a terrible name, hate it. Not elegant
axPAS <- function(nMx, AgeInt, IMR, Sex = "M", region = "W", OAG = TRUE){
	# sex can be "m", "f", or "b"
	# region can be "n","e","s","w",or
	Sex    <- tolower(Sex)
	region <- tolower(region)
	
	N      <- length(nMx)
	ax     <- rep(2.5, N)
	
	if (missing(IMR) | is.na(IMR)){
		IMR <- nMx[1]
	}
	ax[1]  <- geta0CD(M0 = nMx[1], IMR = IMR, Sex = Sex, region = region)
	ax[2]  <- geta1_4CD(M0 = nMx[1], IMR = IMR, Sex = Sex, region = region)
	ax[N]  <- ifelse(!is.na(AgeInt[N]) & AgeInt[N] == 5, 2.5, 1 / nMx[N])
	ax
}

# this is strictly CD west, but what about the others? asked for Preston ref.
# I'd program these differently if I had all the coefs
nMx2ax.greville <- function(nMx, Sex = "m", SRB = 1.05, kludge = FALSE) {
	# QxMx is either Qx or Mx vector of values  
	# with inputQxMx = 1 for Qx as input for QxMx, and 2 for Mx otherwise like in Mortpak LIFETB
	# with Sex ="m" for Male and "f" for Female
	if (Sex == "m"){ # Male under age 5 based on CD-West
		
		#if input is Mx (Table 3.3, p. 48 in Preston et al 2001)
		if (nMx[1] < 0.107) {
			ax0 <- 0.045 + 2.684 * nMx[1]
			ax1 <- 1.651 - 2.816 * nMx[1]
		} else {
			ax0 <- 0.330
			ax1 <- 1.352
		}
		
	}
	if (Sex == "f"){ # Female under age 5 based on CD-West
		#if input is Mx (Table 3.3, p. 48 in Preston et al 2001)
		if (nMx[1] < 0.107) {
			ax0 <- 0.053 + 2.800 * nMx[1]
			ax1 <- 1.522 - 1.518 * nMx[1]
		} else {
			ax0 <- 0.350
			ax1 <- 1.361
		}
	}
	if (Sex == "b"){
		# take SRB weighted avg of the above
	}
	N <- length(nMx)
	## Greville (based on Mortpak LIFETB) for other ages
	
	# this AK jitter replicates a FORTRAN quirk?
	AK             <- 0.1 * log(nMx[3:N] / nMx[1:(N - 2)])
	if (kludge){
	  AK[length(AK)] <- AK[length(AK) -1 ]
    }
	ax <- c(NA, 2.5 - (25 / 12) * (nMx[-c(1, N)] - AK),	1 / nMx[N])
	ax[1:4] <- c(ax0, ax1, 2.5, 2.5)
	if (kludge){
		# Fortran also appears to be doing this, as results cascade...
		ax <- round(ax, 3)
	}
	ax
}

qx2ax.greville <- function(nqx, Sex = "m", SRB = 1.05) {
	# QxMx is either Qx or Mx vector of values  
	# with inputQxMx = 1 for Qx as input for QxMx, and 2 for Mx otherwise like in Mortpak LIFETB
	# with Sex =1 for Male and 2 for Female
	if (Sex=="m"){ # Male under age 5 based on CD-West
		#if input is Qx (Table on p. 20 in Coale and Demeny 1983 for CD-West)
		# TR: need to check these cutpoints, since same value being used for mx and qx...
		if (nqx[1] < 0.1) {
			ax0 <- 0.0425 + 2.875 * nqx[1]
			ax1 <- 1.653 - 3.013 * nqx[1]
		} else 	{
			ax0 <- 0.330
			ax1 <- 1.352
		}
		
	}
	
	if (Sex=="f"){ # Female under age 5 based on CD-West
		#if input is Qx (Table on p. 20 in Coale and Demeny 1983 for CD-West)
		if (nqx[1] < 0.1) {
			ax0 <- 0.050 + 3.000 * nqx[1]
			ax1 <- 1.524 - 1.625 * nqx[1]
		} else 	{
			ax0 <- 0.350
			ax1 <- 1.361
		}
		
	}
	if (Sex == "b"){
		# take SRB weighted avg of the above
	}
	N <- length(nqx)
	## Greville (based on Mortpak LIFETB) for other ages
	ax <- c(NA, 2.5 - (25 / 12) * (nqx[-c(1, N)] - 0.1 * log(nqx[3:N] / nqx[1:(N - 2)])),	
			# TR: this closeout works for mx, but not for qx. Is this just blended out in the iteration?
			1 / nqx[N])
	ax[1:4] <- c(ax0, ax1, 2.5, 2.5)
	ax
}

axUN <- function(nqx, nMx, AgeInt, Sex = "m", tol = .Machine$double.eps, maxit = 1e3, kludge = FALSE){
	stopifnot(!missing(nqx) | !missing(nMx))
	smsq    <- 99999
	
	
	if (missing(nqx) & !missing(nMx)){
# UN (1982) p 31 
# http://www.un.org/esa/population/publications/Model_Life_Tables/Model_Life_Tables.htm
#		For ages 15 and over, the expression for nQx is derived
#		from Greville" as ,nax, = 2.5 - (25/12) (nmx) - k), where
#		k = 1/10 log(nmx+5/nmx-5). For ages 5 and 10, nQx = 2.5
#		and for ages under 5, nQx values from the Coale and
#		Demeny West region relationships are used."
		
		axi <- nMx2ax.greville(nMx = nMx, Sex = Sex, kludge = kludge)
	}
	if (!missing(nqx) & missing(nMx)){
# UN (1982) p 31 
# http://www.un.org/esa/population/publications/Model_Life_Tables/Model_Life_Tables.htm
#		With nqx as input, the procedure is identical, except
#		that an iterative procedure is used to find the nmx and nqx
#		values consistent with the given nqx and with the Greville
#		expression. To complete the life table, the last six nqx
#		values are used to fit the Makeham-type expression
		
		axi <- qx2ax.greville(nqx = nqx, Sex = Sex)
		
		mxi <- mx.from.qx(nqx, ax = axi)
		for (i in 1:maxit) {
			mxi   <- qx2mx(nqx, axi)
			axi   <- Mx2ax.greville(nMx = mxi, Sex = Sex, SRB = SRB)
			qxnew <- mx2qx(mxi, axi)
			smsq  <- sum((qxnew - nqx)^2)
			if (smsq < tol){
				break
			}
		}
	}
	# if both given, then we have ax via identity:
	if (!missing(nqx) & !missing(nMx)){
		axi <- qxmx2ax(nqx = nqx, nmx = nMx)
	}
	
	# if mx, qx, or both are given, then by now we have ax
	axi
}

