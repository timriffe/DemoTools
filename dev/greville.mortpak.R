lt <- readRDS("C:/Users/pgerl/Dropbox/bayesPop/bayesPop2017/LTWorldMale.rds")




greville.mortpak <- function(QxMx, inputQxMx, Sex) {
  # QxMx is either Qx or Mx vector of values  
  # with inputQxMx = 1 for Qx as input for QxMx, and 2 for Mx otherwise like in Mortpak LIFETB
  # with Sex =1 for Male and 2 for Female
 if (Sex==1){ # Male under age 5 based on CD-West
     if (inputQxMx == 1) {  #if input is Qx (Table on p. 20 in Coale and Demeny 1983 for CD-West)
       if (QxMx[1] < 0.1) {
         ax0 <- 0.0425 + 2.875 * QxMx[1]
         ax1 <- 1.653 - 3.013 * QxMx[1]
       } else 
       {
         ax0 <- 0.330
         ax1 <- 1.352
       }
     }
     if (inputQxMx == 2) {  #if input is Mx (Table 3.3, p. 48 in Preston et al 2001)
       if (QxMx[1] < 0.107) {
           ax0 <- 0.045 + 2.684 * QxMx[1]
           ax1 <- 1.651 - 2.816 * QxMx[1]
       } else 
         {
           ax0 <- 0.330
           ax1 <- 1.352
         }
     }
  }
 if (Sex==2){ # Female under age 5 based on CD-West
    if (inputQxMx == 1) {  #if input is Qx (Table on p. 20 in Coale and Demeny 1983 for CD-West)
      if (QxMx[1] < 0.1) {
        ax0 <- 0.050 + 3.000 * QxMx[1]
        ax1 <- 1.524 - 1.625 * QxMx[1]
      } else 
      {
        ax0 <- 0.350
        ax1 <- 1.361
      }
    }
    if (inputQxMx == 2) {  #if input is Mx (Table 3.3, p. 48 in Preston et al 2001)
      if (QxMx[1] < 0.107) {
        ax0 <- 0.053 + 2.800 * QxMx[1]
        ax1 <- 1.522 - 1.518 * QxMx[1]
      } else 
      {
        ax0 <- 0.350
        ax1 <- 1.361
      }
    }
  }
  
  ## Greville (based on Mortpak LIFETB) for other ages
  ax <- c(NA, 2.5 - (25 / 12) * (QxMx[-c(1, length(QxMx))] - 0.1 * log(QxMx[3:length(QxMx)] / QxMx[1:(length(QxMx)-2)])),
            1/QxMx[length(QxMx)])
  ax[1:4] <- c(ax0, ax1, 2.5, 2.5)
  ax
}


greville.mortpak.redux <- function(QxMx, inputQxMx, Sex) {
	
	## Double.MIN_VALUE for log
	DBL_MIN <- .Machine$double.xmin
	
	# n = age interval
	n <- c(1, 4, rep(5, length(QxMx)))
	# with inputQxMx = 1 for Qx as input for QxMx, and 2 for Mx otherwise like in Mortpak LIFETB
	# with Sex =1 for Male and 2 for Female
	if (Sex==1){ # Male under age 5 based on CD-West
		if (inputQxMx == 1) {  #if input is Qx (Table on p. 20 in Coale and Demeny 1983 for CD-West)
			if (QxMx[1] < 0.1) {
				ax0 <- 0.0425 + 2.875 * QxMx[1]
				ax1 <- 1.653 - 3.013 * QxMx[1]
			} else {
				ax0 <- 0.330
				ax1 <- 1.352
			}
		}
		if (inputQxMx == 2) {  #if input is Mx (Table 3.3, p. 48 in Preston et al 2001)
			if (QxMx[1] < 0.1072) {
				ax0 <- 0.045 + 2.684 * QxMx[1]
				ax1 <- 1.651 - 2.816 * QxMx[1]
			} else {
				ax0 <- 0.330
				ax1 <- 1.352
			}
		}
	}
	
	if (Sex==2){ # Female under age 5 based on CD-West
		if (inputQxMx == 1) {  #if input is Qx (Table on p. 20 in Coale and Demeny 1983 for CD-West)
			if (QxMx[1] < 0.1) {
				ax0 <- 0.050 + 3.000 * QxMx[1]
				ax1 <- 1.524 - 1.625 * QxMx[1]
			} else {
				ax0 <- 0.350
				ax1 <- 1.361
			}
		}
		if (inputQxMx == 2) {  #if input is Mx (Table 3.3, p. 48 in Preston et al 2001)
			if (QxMx[1] < 0.1072) {
				ax0 <- 0.053 + 2.800 * QxMx[1]
				ax1 <- 1.522 - 1.518 * QxMx[1]
			} else {
				ax0 <- 0.350
				ax1 <- 1.361
			}
		}
	}
	
	## Initial ax default values
	ax <- c(ax0, ax1, rep(2.5, length(QxMx)-2))
	
	## ax for age 5-9 onward using Greville formula
	for (j in 3:(length(QxMx)-1)) {
		## Mortpak LIFETB for age 5-9 and 10-14
		# ax[j] = 2.5 
		## for ages 15-19 onward
		## AK <- log(QxMx[j+1]/QxMx[j-1])/10
		## ax[j] <- 2.5 - (25.0/12.0) * (QxMx[j] - AK)
		
		## improved Greville formula for adolescent ages 5-9 and 10-14
		## Let the three successive lengths be n1, n2 and n3, the formula for 5a5 is:
		## ax[i] = 2.5 - (25 / 12) * (mx[i] - log(mx[i + 1] / mx[i-1])/(n1/2+n2+n3/2))
		## for age 5-9, coefficient should be 1/9.5, because age group 1-4 has only 4 ages (not 5), while the other 5-year age group are 1/10
		## ax[i] = 2.5 - (25 / 12) * (mx[i] - (1/9.5)* log(mx[i + 1] / mx[i-1]))
		## Age 20-25, ..., 95-99
		## Greville (based on Mortpak LIFETB) for other ages, new implementation
		ax[j] = (n[j]/2) - ((n[j]^2) / 12) * (QxMx[j] - log(max(QxMx[j + 1] / max(QxMx[j - 1], DBL_MIN), DBL_MIN)) / (n[j-1]/2 + n[j] + n[j+1]/2))
		
		
		## add constraint at older ages (in Mortpak and bayesPop)
		## 0.97 = 1-5*exp(-5)/(1-exp(-5)), for constant mu=1, Kannisto assumption  (Mortpak used 1 instead of 0.97).
		if (j > 9 && ax[j] < 0.97) {ax[j] = 0.97}
		
		## Extra condition based on Mortpak LIFETB for age 65 onward
		if (inputQxMx == 1 & j >= 15) {
			tmp <- 0.8 * ax[j-1]
			if (ax[j] < tmp) {ax[j] <- tmp}
		}
		
	}
	ax[length(QxMx)] <- 1/QxMx[length(QxMx)]
	if (ax[length(QxMx)] < 0.97) {ax[length(QxMx)] = 0.97}
	if (inputQxMx == 1) {
		tmp <- 0.8 * ax[length(QxMx)-1]
		if (ax[length(QxMx)] < tmp) {ax[length(QxMx)]  <- tmp}
	}
	
	return(ax)
}

QxMx <- qx



mx.from.qx <- function(qx, ax, n = inferAgeIntAbr(vec = qx)) {
    qx / (n - (n - ax) * qx)      }

qx.from.mx <- function(mx, ax, n = inferAgeIntAbr(vec = qx)) {
    (n * mx) / (1.0 + ((n - ax) * mx)) 
}

ax.from.qxmx <- function(qx, mx, n = inferAgeIntAbr(vec = qx)){
	1 / mx - n / qx + n
}

qx <- lt$qx

## initialize ax using Qx
ax1 <- greville.mortpak(qx, 1, 1)
axi <- ax1

## iterate until output qx match input qx using mxi with axi and qx as input
smsq <- 99999
epsilon <- 1.00E-16
while(smsq > epsilon) {
  mxi   <- mx.from.qx(qx, axi, n)
  axi   <- greville.mortpak(mxi, 2)
  qxnew <- qx.from.mx(mxi, axi, n)
  smsq  <- sum((qxnew - qx)^2)
  print(smsq)
}

mx <- mxi
ax <- axi
## note the last mx for the open age group needs to be computed differently, 
# or rather mx only up to age 95-99 used to fit Kannisto, and the open age group discarded.
mx[length(mx)] <- NA
ax[length(mx)] <- NA

check <- data.frame(cbind(mx, qxnew, qx, ax))

#install.packages("bayesPop")

# this is strictly CD west, but what about the others? asked for Preston ref.
# I'd program these differently if I had all the coefs
Mx2ax.greville <- function(nMx, Sex = "m", SRB = 1.05) {
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
	ax <- c(NA, 2.5 - (25 / 12) * (nMx[-c(1, N)] - 0.1 * log(nMx[3:N] / nMx[1:(N - 2)])),	1 / nMx[N])
	ax[1:4] <- c(ax0, ax1, 2.5, 2.5)
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



axUN <- function(qx, mx, Sex, tol = .Machine$double.eps, maxit = 1e3){
	stopifnot(!missing(qx) | !missing(mx))
	smsq    <- 99999
	

	if (missing(qx) & !missing(mx)){
# UN (1982) p 31 
# http://www.un.org/esa/population/publications/Model_Life_Tables/Model_Life_Tables.htm
#		For ages 15 and over, the expression for nQx is derived
#		from Greville" as ,nax, = 2.5 - (25/12) (nmx) - k), where
#		k = 1/10 log(nmx+5/nmx-5). For ages 5 and 10, nQx = 2.5
#		and for ages under 5, nQx values from the Coale and
#		Demeny West region relationships are used."
		
		axi <- mx2ax.greville(nMx = mx, Sex = Sex)
	}
	if (!missing(qx) & missing(mx)){
# UN (1982) p 31 
# http://www.un.org/esa/population/publications/Model_Life_Tables/Model_Life_Tables.htm
#		With nqx as input, the procedure is identical, except
#		that an iterative procedure is used to find the nmx and nqx
#		values consistent with the given nqx and with the Greville
#		expression. To complete the life table, the last six nqx
#		values are used to fit the Makeham-type expression
		
		axi <- qx2ax.greville(nqx = qx, Sex = Sex)
		mxi <- mx.from.qx(qx, ax = axi)
		for (i in 1:maxit) {
			mxi   <- mx.from.qx(qx, axi)
			axi   <- Mx2ax.greville(mxi, Sex = Sex, SRB = SRB)
			qxnew <- qx.from.mx(mxi, axi)
			smsq  <- sum((qxnew - qx)^2)
			if (smsq < tol){
				break
			}
		}
	}
	# if both given, then we have ax via identity:
	if (!missing(qx) & !missing(mx)){
		axi <- ax.from.qxmx(qx = qx, mx = mx)
	}
	
	# if mx, qx, or both are given, then by now we have ax
	axi
}



