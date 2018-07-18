
#TODO in top to CENSUR.R: check that data in 5 year ages, group to 5-year ages if not before running.
# Author: tim
###############################################################################
#' Census survival estimation
#' 
#' @description This function family implements the spreadsheets from the UN repository 
#' "Adult Mortality Estimation Files/countries/JAPAN/CENSUS~1.SUR", prepared by Griffith Feeney ca 2007.
#' 
#' @param pop1 numeric. Vector of population counts in 5-year age groups of census 1.
#' @param pop2 numeric. Vector of population counts in 5-year age groups of census 2.
#' @param Age numeric. Vector of ages corresponding to the lower integer bound of the counts.
#' @param date1 date. Date of the first census. See details for ways to express it.
#' @param date2 date. Date of the second census. See details for ways to express it.
#' @param exprior numeric. Vector of remaining life expectancy (not e(0), see details) for same population and time but from a different source.

#' @details The lengths of \code{pop1} and \code{pop2} should match. We assume that the last age group is open, and throw it out. 
#' If your last age group is not open and you want to keep it, append an element, such as \code{NA} to the end of your vectors. 
#' Dates can be given in three ways 1) a \code{Date} class object, 2) an unambiguous character string in the format \code{"YYYY-MM-DD"}, 
#' or 3) as a decimal date consisting in the year plus the fraction of the year passed as of the given date. Value for life expectancy at birth e(0) in \code{exprior} 
#' should be set to \code{NA}.
#' @references 
#' \insertRef{united1955manual}{DemoTools}
#' \insertRef{united1983manual}{DemoTools}
#' @return The median absolute percent error of the census survival remaining life expectancy vector compared with the external remaining life expectancy estimate (\code{exprior}).
#' @export
#' @importFrom stats median
#' @examples
#'# Based on Feeney's spreadsheets CENSUR~1.XLS, CENSUR~2.XLS, and CENSUR~3.XLS
#'# this is census data from Japan
#' 
#' # create a bunch of data to show the three variants.
#' pop1960 <- c(3831870,4502304,5397061,4630775,4193184,4114704,3770907,3274822,
#'    2744786,2559755,2160716,1839025,1494043,1133409,870238,577972,313781,131547)
#' pop1965 <- c(3983902,3854281,4513237,5373547,4572392,4206801,4110076,3751030,
#'    3231736,2697217,2485095,2071540,1719370,1343444,955567,644043,341170,176068)
#' pop1970 <- c(4292503,3988292,3852101,4492096,5347327,4571868,4190340,4085338,
#'    3674127,3198934,2648360,2382691,1970485,1584699,1172155,736258,408191,53116)
#' Age          <- seq(0, 85, by = 5)
#' date1960     <- "1960-10-01"
#' date1965     <- "1965-10-01"
#' date1970     <- "1970-10-01"
#' date1960fake <- "1960-10-05"
#' ex1  <- c(NA,69.45,64.62,59.72,54.87,50.11,45.37,40.65,35.99,
#' 31.40,26.94,22.64,18.52,14.67,11.20,8.25,5.95,6.00)
#' ex2  <- c(NA,70.19,65.33,60.41,55.54,50.74,45.96,41.21,
#' 36.52,31.89,27.39,23.05,18.89,14.99,11.45,8.43,6.06)
#' 
#' # reproduces CENSUR~1.XLS
#' # this method called under the hood when censuses are a clean 10 years apart
#' survRatioError(
#'  pop1 = pop1960, 
#'  pop2 = pop1970, 
#'  Age = Age, 
#'  date1 = date1960, 
#'  date2 = date1970, 
#'  exprior = ex1) #.44
#' 
#' # reproduces CENSUR~3.XLS
#' # this method called under the hood when censuses are a clean 5 years apart
#' survRatioError(
#' 		pop1 = pop1965,
#' 		pop2 = pop1970,
#' 		Age = Age,
#' 		date1 = date1965,
#' 		date2 = date1970,
#' 		exprior = ex2) #.40
#' 
#' # reproduces CENSUR~2.XLS (if indeed censuses are oddly spaced)
#' # this method called under the hood when censuses are not spaced 
#' # a clean 5 or 10 years apart.
#' survRatioError(
#' 		pop1 = pop1960,
#' 		pop2 = pop1970,
#' 		Age = Age,
#' 		date1 = date1960fake,
#' 		date2 = date1970,
#' 		exprior = ex1)# 1.067818

# result matches CENSUR~2 exactly if you change the
# decimal dates used in the spreadsheet to the following
# print(dec.date(date1960fake),digits=15)
# can't call up the synthetic method on the original dates
# because it gets automatically directed to the cleaner 10-year
# staggered method, which is preferred.

survRatioError <- function(pop1, pop2, Age, date1, date2, exprior){
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
	N       <- N - 1
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
	
	# Lx approximated as avg of lx
	Lx         <- c(NA, ma(lx, 2))
	# 2.5 times the pairwise sums of that to get the 5-year Lx approximation
	Lx5        <- c(2.5 * (Lx[-length(Lx)] + Lx[-1]))
	Lx5[N]     <- Lx[N] * exprior[N]
	
	# cut to size
	Lx         <- Lx[1:N]
	Lx5        <- Lx5[1:N]
	exprior    <- exprior[1:N]
	
	Tx         <- rev(cumsum(rev(Lx5)))
	ex.est     <- Tx / Lx
	dif        <- exprior - ex.est
	abserror   <- 100 * abs(dif / exprior)
	abserror[N] <- median(abserror[-N], na.rm = TRUE)
	
	# return the median absolute percent error
	median(abserror, na.rm = TRUE)
}

#' Estimate survival curve from censuses spaced 10 years apart.
#' 
#' @description This function reproduces calculations through column H of \code{CENSUR~1.XLS} by 
#' Griff Feeney. Censuses are assumed to be spaced 10 years apart and population counts for both censuses in 5-year age groups. 
#' The staggered ratio of these approximates lifetable function npx =1-nqx, i.e. probability of surviving between age x and age x+n. 
#' The cumulative product of this approximates the survival function lx.  
#' 
#' @param pop1 numeric. Vector of population counts in 5-year age groups of census 1.
#' @param pop2 numeric. Vector of population counts in 5-year age groups of census 2.
#' @return Survival function lx as a numeric vector with radix 1.
#' @details Checking for census spacing of 10 years must happen prior to this function. 
#' Also, we assume the open age group has already been trimmed off.
#' @export

#' @examples
#'# 1960 vs 1970 pops
#' pop1 <- c(3831870,4502304,5397061,4630775,4193184,4114704,3770907,3274822,
#' 		2744786,2559755,2160716,1839025,1494043,1133409,870238,577972,313781)
#' pop2 <- c(4292503,3988292,3852101,4492096,5347327,4571868,4190340,4085338,
#' 		3674127,3198934,2648360,2382691,1970485,1584699,1172155,736258,408191)
#' surv10(pop1,pop2)

surv10 <- function(pop1, pop2){
	N       <- length(pop1)
	stagger <- 2
	ratio   <- pop2[-c(1:stagger)] / pop1[1:(N - stagger)]
	
	NN       <- length(ratio)
	sx1     <- 1:NN %% 2 == 1
	sx2     <- 1:NN %% 2 == 0
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
	lx[1:N]
}

#' Estimate survival curve from censuses spaced 5 years apart.
#' 
#' @description This function reproduces calculations of \code{CENSUR~3.XLS} by 
#' Griff Feeney. Censuses are assumed to be spaced 5 years apart and population counts for both censuses in 5-year age groups. 
#' The staggered ratio of these approximates lifetable function npx =1-nqx, i.e. probability of surviving between age x and age x+n. 
#' The cumulative product of this approximates the survival function lx.  
#' 
#' @param pop1 numeric. Vector of population counts in 5-year age groups of census 1.
#' @param pop2 numeric. Vector of population counts in 5-year age groups of census 2.
#' @return Survival function lx as a numeric vector with radix 1.
#' @details Checking for census spacing of 5 years must happen prior to this function. 
#' Also, we assume the open age group has already been trimmed off.
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

#' Estimate survival curve from censuses spaced N years apart.
#' 
#' @description This function reproduces calculations of \code{CENSUR~2.XLS} by 
#' Griff Feeney. Censuses are assumed to be spaced N years apart and population counts for both censuses in 5-year age groups. 
#' The staggered ratio of these approximates lifetable function npx =1-nqx, i.e. probability of surviving between age x and age x+n. 
#' The cumulative product of this approximates the survival function lx.  
#' 
#' @param pop1 numeric. Vector of population counts in 5-year age groups of census 1.
#' @param pop2 numeric. Vector of population counts in 5-year age groups of census 2.
#' @param interval numeric. Length of the inter-censal period.
#' 
#' @return Survival function lx as a numeric vector with radix 1.
#' @details Checking for census spacing of N years must happen prior to this function. 
#' Also, we assume the open age group has already been trimmed off.
#' @export
#' @details We assume the open age group has already been trimmed off. 
#' The time interval N must be measured as a decimal in advance
#' and specified. This method uses a synthetic approximation of person-years lived in each age interval
#' over the intercensal period and then a second approximation based on age-specific growth rates to 
#' produce an estimate of lifetable px. This value of px is not bounded to [0,1], and therefore the 
#' resulting lx approximation is not strictly positive or monotonically non-increasing, 
#' so downstream usage of this result is limited. For example, is not advisable to use it in the 
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

