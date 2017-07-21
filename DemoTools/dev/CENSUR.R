
# Author: tim
###############################################################################

# based on Feeney's spreadsheet CENSUR~1.XLS
pop1 <- c(3831870,4502304,5397061,4630775,4193184,4114704,3770907,3274822,
2744786,2559755,2160716,1839025,1494043,1133409,870238,577972,313781,131547)
pop2 <- c(4292503,3988292,3852101,4492096,5347327,4571868,4190340,4085338,
3674127,3198934,2648360,2382691,1970485,1584699,1172155,736258,408191,53116)
Age  <- seq(0,85,by=5)
date1 <- "1960-10-01"
date2 <- "1970-10-01"

survRatio <- function(pop1, pop2, Age, date1, date2, exprior){
	stopifnot(all.equal(length(pop1),length(pop2),length(Age)))
	
	# census spacing used to stagger ages.
	date1   <- dec.date(date1)
	date2   <- dec.date(date2)
	# calculate intercensal period
	int     <- date2 - date1
	
	# assume final age interval is open, disregard it
	N       <- length(pop1)
	pop1    <- pop1[-N]
	pop2    <- pop2[-N]
	Age     <- Age[-N]
	
	# this 
	if (abs(int - 5) < .01){
		# are we super close to 5 years apart?
		lx <- surv5(pop1, pop2)
	} else {
	  if (abs(int - 10) < .01){
		# are we super close to 10 years apart?
		lx <- surv10(pop1, pop2)
	  } else {
		# otherwise use the synthetic method
		lx <- survN(pop1, pop2, int)
	  }
    }
	# all remaining columns should be identical.

	
}

#' estimate survival curve from censuses spaced 5 years apart
#' 
#' @description This simple function reproduces calculations through column H of \code{CENSUR~1.XLS} by 
#' Griff Feeney. We assume censuses spaced 10 years apart and population counts for both censuses in 5-year age groups. 
#' The staggered ratio of these is like the lifetable px (1-qx). The cumulative product of this then spits back something 
#' we can use to approximate lx.  
#' 
#' @param pop1 numeric vector population counts of census 1
#' @param pop2 numeric vector population counts of census 2
#' @return lx numeric vector of surivorship (radix 1)
#' @details Checking for census spacing of 10 years must happen prrior to this function. 
#' Also, we assume the open age group has already been trimmed off. No checking done here at all.
#' @export

#' @examples
#'# 1960 vs 1970 pops
#' pop1 <- c(3831870,4502304,5397061,4630775,4193184,4114704,3770907,3274822,
#' 		2744786,2559755,2160716,1839025,1494043,1133409,870238,577972,313781)
#' pop2 <- c(4292503,3988292,3852101,4492096,5347327,4571868,4190340,4085338,
#' 		3674127,3198934,2648360,2382691,1970485,1584699,1172155,736258,408191)
#' surv10(pop1,pop2)

surv10 <- function(pop1, pop2){
	stagger <- 2
	ratio   <- pop2[-c(1:stagger)] / pop1[1:(N-stagger-1)]
	
	N       <- length(ratio)
	sx1     <- 1:N %% 2 == 1
	sx2     <- 1:N %% 2 == 0
	lx1     <- c(1, cumprod(ratio[sx1]))
	lx2     <- mean(lx1[1:2]) * c(1, cumprod(ratio[sx2]))
	
	# at this point lx vectors can be of different lengths.
	# need to be of same length for interleaving trick to work.
	N1 <- length(lx1)
	N2 <- length(lx2)
	if (N2 < N1){
		lx2 <- c(lx2, rep(NA, N1 - N2))
	}
	# use object trick to interleave
	lx      <- c(rbind(lx1, lx2))
	lx
}

#' estimate survival curve from censuses spaced 5 years apart
#' 
#' @description This simple function reproduces calculations through column F of \code{CENSUR~3.XLS} 
#' by Griff Feeney. We assume censuses spaced 5 years apart and population counts for both censuses in 5-year age groups. 
#' The staggered ratio of these is like the lifetable px (1-qx). The cumulative product of this then spits back something 
#' we can use to approximate lx.  
#' 
#' @param pop1 numeric vector population counts of census 1
#' @param pop2 numeric vector population counts of census 2
#' @return lx numeric vector of surivorship (radix 1)
#' @details Checking for census spacing of 5 years must happen prrior to this function. 
#' Also, we assume the open age group has already been trimmed off. No checking done here at all.
#' @export

#' @examples
#'# 1965 vs 1970 pops
#' pop1 <- c(3983902,3854281,4513237,5373547,4572392,4206801,4110076,3751030,
#' 3231736,2697217,2485095,2071540,1719370,1343444,955567,644043,341170)
#' pop2 <- c(4292503,3988292,3852101,4492096,5347327,4571868,4190340,4085338,
#' 3674127,3198934,2648360,2382691,1970485,1584699,1172155,736258,408191)
#' surv5(pop1,pop2)
surv5 <- function(pop1, pop2){
	N       <- length(pop1)
	stagger <- 1
	ratio   <- pop2[-c(1:stagger)] / pop1[1:(N - stagger)]
	lx      <- cumprod(c(1,ratio))
	lx
}

#' estimate survival curve from censuses spaced N years apart
#' 
#' @description This simple function reproduces calculations through column H of \code{CENSUR~2.XLS} 
#' by Griff Feeney. We have censuses spaced an arbitrary N years apart and population counts for 
#' both censuses in 5-year age groups. 
#' @param pop1 numeric vector population counts of census 1
#' @param pop2 numeric vector population counts of census 2
#' @param interval a numeric value, annualized intercensal period.
#' @return lx numeric vector of surivorship (radix 1)
#' @details We assume the open age group has already been trimmed off. No checking done here at all. 
#' The time interval N must be measured as a decimal in advance
#' and specified. This method uses a synthetic approximation of person-years lived in each age interval
#' over the intercensal period and then a second approximation based on age-specific growth rates to 
#' produce an estimate of lifetable px. This value of px is not bounded to [0,1], and therefore the 
#' resulting lx approximation is not strictly positive or monotonically non-increasing, 
#' so downstream usage of this result is limited. For example, you wouldn't want to use it in the 
#' standard way to derive the lifetable dx.
#' 
#' @export

#' @examples
#'# 1960 vs 1970 pops
#' pop1 <- c(3831870,4502304,5397061,4630775,4193184,4114704,3770907,3274822,
#' 		2744786,2559755,2160716,1839025,1494043,1133409,870238,577972,313781)
#' pop2 <- c(4292503,3988292,3852101,4492096,5347327,4571868,4190340,4085338,
#' 		3674127,3198934,2648360,2382691,1970485,1584699,1172155,736258,408191)
#' interval <- 10
#' survN(pop1,pop2,interval)

survN <- function(pop1, pop2, interval){
	r   <- log(pop2 / pop1) / interval
	PYL <- (pop2 - pop1) / (interval * r)
	Sx  <- (PYL[-1] * exp(2.5*r[-1])) / (PYL[-length(PYL)] * exp(-2.5 * r[-length(r)]))
	lx  <- cumprod(c(1, Sx))
	lx
}
