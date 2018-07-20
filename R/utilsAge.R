
# Author: tim
###############################################################################
# contains utilities having to do with age groups.
# - simple operations like detecting, grouping, and so forth.
# not complex things like ungrouping ages.

#' trim two age vectors to matching N-year age groups
#' 
#' @description determining which N-year (e.g. 5-year) age groups two vectors have in common
#' is helpful for matching vector lengths, and questions of alignment. Used as a utility throughout.

#' @param Age1 integer vector of first age groups (lower bounds)
#' @param Age2 integer vector of second age groups (lower bounds)
#' @param N integer target age group interval (just one number)
#' @param consecutive logical default \code{TRUE}. Throw error if resulting age groups not consecutive?
#' @param ageMin integer optional lower age bound for output
#' @param ageMax integer optional upper age bound for output
#' 
#' @return an integer vector of the N-year age groups present in both \code{Age1} and \code{Age2}
#' 
#' @export 
#' 
#' @examples 
#' Age1 <- seq(0, 100, by = 5)
#' Age2 <- 0:80
#' AGEN(Age1, Age2, N = 5)
#' AGEN(Age1, Age2, N = 5, ageMax = 70)
AGEN <- function(Age1, Age2, N = 5, consecutive  = TRUE, ageMin = 0, ageMax = max(c(Age1,Age2))){
	age1_5 <- Age1[Age1 %% 5 == 0]
	age2_5 <- Age2[Age2 %% 5 == 0]
	
	# ages in common only
	ageN <- sort(intersect(age1_5, age2_5))
	ageN <- ageN[ageN >= ageMin & ageN <= ageMax]
	
	# make sure it's consecutive series
	if (consecutive){
		stopifnot(all(diff(ageN) == N))
	}
	
	ageN
}

#' calculate which large age group single ages belong to
#' 
#' @description Assign single ages to age groups of equal and arbitrary width, and also optionally shifted.
#' 
#' @param Age integer vector of single ages (lower bound)
#' @param N integer. Desired width of resulting age groups
#' @param shiftdown integer (optional), move the grouping down by one or more single ages?
#' 
#' @details If you shift the groupings, then the first age groups may have a negative lower bound 
#' (for example of -5). These counts would be discarded for the osculatory version of Sprague smoothing,
#' for example, but they are preserved in this function. The important thing to know is that if you shift 
#'  the groups, the first and last groups won't be N years wide. For example if \code{shiftdown} is 1, 
#' the first age group is 4-ages wide. 
#'  
#' @return an integer vector of \code{length(Age)} indicating the age group that each single age belongs to.
#' 
#' @export
#' @examples
#' Age <- 0:100
#' calcAgeN(Age)
#' calcAgeN(Age, N = 4)
#' calcAgeN(Age, N = 3)
#' calcAgeN(Age, shiftdown = 1)
#' calcAgeN(Age, shiftdown = 2)
#' # can be used to group abridged into equal 5-year intervals
#' AgeAbr <- c(0,1,5,10,15,20)
#' calcAgeN(AgeAbr)
calcAgeN <- function(Age, N = 5, shiftdown = 0){
	shift <- abs(shiftdown)
	stopifnot(shift < N)
	(Age + shift) - (Age + shift) %% N 
}

#' calculate which abridged age group single ages belong to
#' 
#' @description Assign single ages to 5-year abridged age groups. That means that age 0 is kept as a single age,
#' ages 1-4 are grouped together as abridged age 1, and thereafter 5-year age groups are used.
#' 
#' @param Age integer vector of single ages (lower bound)
#' 
#' @details In the case that the single \code{Age} vector starts at some age higher than 4, 
#' this is just the same as \code{calcAgeN(Age,5,0)}.
#'  
#' @return an integer vector of \code{length(Age)} indicating the abridged age group that each single age belongs to.
#' 
#' @export
#' @examples
#' Age <- 0:70
#' calcAgeAbr(Age)
#' calcAgeN(Age,5,0)
calcAgeAbr <- function(Age){
	Abr               <- Age - Age %% 5
	Abr[Age %in% 1:4] <- 1
	Abr
}

#' infer abridged age groups widths
#' 
#' @description This function is an auxiliary used by top level functions where it is 
#' guaranteed that age groups are standard abridged age groups. If \code{Age} is specified,
#' this will work as well for other age groupings.
#' 
#' @param Age integer vector. Lower bound of each age group.
#' @param vec any vector, presumably a count, rate, or similar.
#' @param OAG logical (default \code{FALSE}). Is the final age group open?
#' @param OAGvalue numeric or integer. The value to use for the final age interval if \code{OAG = TRUE}. Default \code{NA}.
#' 
#' @details If based soley on the length of a vector, this will give erroneous results if ages 
#' are anything other than standard abridged ages groups. If the final age group is open, the 
#' interval with is defined as \code{NA}. \code{Inf} or \code{-1} would have 
#' also been a good choice, but we went with \code{NA}.
#' 
#'  
#' @return an integer vector of \code{length(vec)} indicating the width of the abridged age group that each 
#' vector element corresponds to.
#' 
#' @export
#' @examples
#' vec <- runif(20)
#' inferAgeIntAbr(vec = vec)
#' inferAgeIntAbr(vec = vec, OAG = TRUE)
inferAgeIntAbr <- function(Age, vec, OAG = FALSE, OAGvalue = NA){
	
	# Age is preferred (lower bounds)
	if (!missing(Age)){
		ageint <- age2int(Age = Age, 
				          OAG = OAG, 
						  OAGvalue = OAGvalue)
	}
	# otherwise length of vector will do
	if (missing(Age) & !missing(vec)){
		stopifnot(length(vec)>5)
		ageint    <- rep(5,length(vec))
		ageint[1] <- 1
		ageint[2] <- 4
		if (OAG){
			ageint[length(ageint)] <- OAGvalue
		}
	}
	
	ageint
}

#' determine abridged ages up to a given maximum age group
#' @description Produce standard abridged age groups (lower bounds) up to a specified maximum age group. 
#' @details If the highest age group is not evenly divisible by 5 then age classes only go up to its 5-year lower bound.
#' @param maxA integer. Default 80
#' @return integer vector of ages \code{c(0,1,5,10,15,...)}
#' @export 
#' @examples 
#' maxA2abridged(80)
#' stopifnot(all(maxA2abridged(100) == maxA2abridged(102)))
maxA2abridged <- function(maxA = 80){
	sort(unique(calcAgeAbr(0:maxA)))
}



#' infer lower age bounds from age class intervals
#' @description Determine lower bounds of age classes based on a vector of age intervals and a starting age.
#' 
#' @param AgeInt integer or numeric vector of age intervals 
#' @param minA the lowest age, default 0.
#' @export 
#' @return Age vector, of same length as \code{AgeInt}.
#' @examples 
#' AgeInt <- c(1,4,rep(5,17))
#' int2age(AgeInt)
int2age <- function(AgeInt, minA = 0){
	cumsum(AgeInt) - AgeInt + minA
}

#' infer age class intervals from lower age bounds
#' @description Determine age calss intervals based on a vector of age class lower bounds.
#' @details If the final age group is open, it is given a value of \code{NA} by default, or else a user-determined value. 
#' If the final age group is closed, it is assumed to be equal to the next-lower interval. If the final age interval is 
#' known and not equal to the next lowest interval, specify \code{OAG = TRUE} and assign its value to \code{OAvalue}.
#' @param Age integer or numeric vector of age intervals 
#' @param OAG logical. \code{TRUE} if the final age group is open.
#' @param OAGvalue numeric or integer. The value to use for the final age interval if \code{OAG = TRUE}. Default \code{NA}.
#' @export 
#' @return Age interval vector, of same length as \code{Age}.
#' @examples 
#' # single age examples:
#' Age  <- 0:100
#' age2int(Age, OAG = TRUE, OAvalue = NA)
#' age2int(Age, OAG = TRUE, OAvalue = Inf)
#' age2int(Age, OAG = FALSE)
#' 
#' # and for abridged ages:
#' AgeA <- c(0,1,seq(5,85,by=5))
#' age2int(AgeA, OAG = TRUE, OAvalue = NA)
#' age2int(AgeA, OAG = TRUE, OAvalue = Inf)
#' age2int(AgeA, OAG = FALSE)
age2int <- function(Age, OAG = TRUE, OAvalue = NA){
	fd <- diff(Age)
	c(fd,ifelse(OAG,OAvalue,fd[length(fd)]))
}


#' group single ages into equal age groups of arbitrary width
#' 
#' @description This can be useful to check constrained sums, or as an intermediate step for smoothing.
#' 
#' @param Value numeric vector of single age counts
#' @param Age integer vector of lower bounds of single age groups
#' @param N integer. Default 5. The desired width of resulting age groups
#' @param shiftdown integer. Default 0. Optionally shift age groupings down by single ages 
#' @param AgeN optional integer vector, otherwise calculated using \code{calcAgeN()}
#' @param OAnew integer value of lower bound of new open age group
#' @return vector of counts in N-year age groups
#' 
#' @details If you shift the groupings, then the first age groups may have a negative lower bound 
#' (for example of -5). These counts would be discarded for the osculatory version of Sprague smoothing,
#' for example, but they are preserved in this function. The important thing to know is that if you shift 
#'  the groups, the first and last groups won't be N years wide. For example if \code{shiftdown} is 1, the first age group is 4-ages wide. The ages themselves are not returned, 
#' but they are the name attribute of the output count vector. Note this will also correctly group abridged ages
#' into equal 5-year age groups if the \code{Age} argument is explicitly given. \code{OAnew} (optional) must be less than or equal to \code{max(Age)} to have any effect.
#' @export
#' @examples
#'  India1991males <- c(9544406,7471790,11590109,11881844,11872503,12968350,11993151,10033918
#'  		,14312222,8111523,15311047,6861510,13305117,7454575,9015381,10325432,
#'           9055588,5519173,12546779,4784102,13365429,4630254,9595545,4727963
#'  		,5195032,15061479,5467392,4011539,8033850,1972327,17396266,1647397,
#'           6539557,2233521,2101024,16768198,3211834,1923169,4472854,1182245
#'  		,15874081,1017752,3673865,1247304,1029243,12619050,1499847,1250321,
#'           2862148,723195,12396632,733501,2186678,777379,810700,7298270
#'  		,1116032,650402,1465209,411834,9478824,429296,1190060,446290,362767,
#'           4998209,388753,334629,593906,178133,4560342,179460
#'  		,481230,159087,155831,1606147,166763,93569,182238,53567,1715697,
#'           127486,150782,52332,48664,456387,46978,34448
#'  		,44015,19172,329149,48004,28574,9200,7003,75195,13140,5889,18915,21221,72373)
#'  Age <- 0:100
#'  groupAges(India1991males, N = 5)
#'  groupAges(India1991males, N = 5, shiftdown = 1)
#'  groupAges(India1991males, N = 5, shiftdown = 2)
#'  groupAges(India1991males, N = 5, shiftdown = 3)
#'  groupAges(India1991males, N = 5, shiftdown = 4)
#'  groupAges(India1991males, N = 5, OAnew = 80) 
 
groupAges <- function(Value, 
		Age = 1:length(Value) - 1, 
		N = 5, 
		shiftdown = 0, 
		AgeN,
		OAnew = max(Age)){
	if (missing(AgeN)){
		AgeN <- calcAgeN(Age = Age, N = N, shiftdown = shiftdown)
	}
	out <- tapply(Value, AgeN, sum)
	
	# group down open age 
	if (OAnew < max(AgeN)){
		AgeA <- sort(unique(AgeN))
		out  <- groupOAG(Value = out, Age = AgeA, OAnew = OAnew)
	}
	out
}

#' Group down to a new open age class
#' 
#' @description This simple utility lowers the open age group. It only returns the input value vector, not the age vector.
#' @param Value numeric vector of counts
#' @param Age integer vector of age classes
#' @param OAnew the desired open age group
#' @export 
#' @return Value vector potentially of reduced length up to OAG

groupOAG <- function(Value, Age, OAnew){
	stopifnot(OAnew %in% Age)
	N        <- length(Value[Age <= OAnew])
	Value[N] <- sum(Value[Age >= OAnew])
	Value    <- Value[1:N]
	Value
}


#' logical checking of whether age classes appear single
#' 
#' @description check whether a vector of ages consists in single ages. This
#' makes sense sometimes when age intervals are not explicitly given.
#' 
#' @param Age integer vector of age classes
#' 
#' @return logical \code{TRUE} if detected as single ages, \code{FALSE} otherwise
#' 
#' @details In cases where ages are indeed single, but repeated, this will still return \code{FALSE}. 
#' Therefore make sure that the age vector given refers to a single year of a single population, etc.
#' @export
#' @examples
#' Age <- 0:99
#' Age2 <- c(0:10,0:10)
#' Age3 <- seq(0,80,by=5)
#' Age4 <- seq(0,10,by=.5)
#' is.single(Age)  # TRUE
#' is.single(Age2) # FALSE repeated, can't tell.
#' is.single(Age3) # FALSE not single ages
#' is.single(Age4) # FALSE not single ages
is.single <- function(Age){
	all(diff(sort(Age)) == 1)
}

#' detect if a vector of lower age bounds is plausibly of abridged ages
#' @description a logical utility that checks if a given vector is of the lower bounds of abridged age groups or not.
#' @param Age integer vector of lower age bounds
#' @export
#' @return logical \code{TRUE} if abridged, \code{FALSE} otherwise.

#' @examples
#' # as expected, TRUE
#' is.abridged(c(0,1,5,10,15,20,25))
#' # standard 5, not abridged, FALSE
#' is.abridged(c(0,5,10,15,20,25))
#' # plausible, TRUE
#' is.abridged(c(1,5,10,15,20,25))
#' # plausible, TRUE
#' is.abridged(c(5,10,15,20,25))
#' # 10 year age group not abridged, FALSE
#' is.abridged(c(0,1,5,10,15,25))
is.abridged <- function(Age){
	maxA        <- max(Age)
	minA        <- min(Age)
	abr_default <- calc_AgeAbr_OAG(maxA)
	abr_default <- abr_default[abr_default >= minA]
	identical(Age,abr_default)
}

# deprecated functions

# TR: deprecated July 20, 2018
##' infer lower age class bounds from vector of age intervals
##' 
##' @description a simple identity
##' @details It makes no difference whether the final age group is open or how that is coded,
##' as long as all lower age intervals are coercible to integer.
##' @param AgeInt vector of age intervals
##' @param minA what is the lowest age bound to start counting from (default 0)?
##' @return integer vector of ages
##' @export 
##' @examples 
##' int1 <- c(rep(1,100),NA)
##' int2 <- c(rep(1,100),"+")
##' int3 <- inferAgeIntAbr(vec = rep(1,20), OAG = TRUE)
##' AgefromAgeInt(int1)
##' AgefromAgeInt(int2) # character OK
##' AgefromAgeInt(int3)
## TR: does this give upper bounds rather than lower bounds
## of age?
#AgefromAgeInt <- function(AgeInt, minA = 0){
#	N      <- length(AgeInt)
#	AgeInt <- as.integer(AgeInt[-N])
#	cumsum(c(minA, AgeInt))
#} 

# TR deprecated 20 July, 2018. Use groupAges() instead

##' aggregates single year age groups into 5 year age groups
##' 
##' @description Creates five year age groups from single year ages. 
##' @details Sums five year age intervals
##' 
##' @param Value numeric vector of single year age groups.
##' 
##' @export 
##' @examples 
##' MalePop <- seq(1,100)
##' convertSingleTo5Year(MalePop)
#
#convertSingleTo5Year <- function(Value){
#	shiftZero  <- Value
#	shiftOne   <- Value[-1]
#	shiftTwo   <- shiftOne[-1]
#	shiftThree <- shiftTwo[-1]
#	shiftFour  <- shiftThree[-1]
#	
#	# TR: not sure what to make of zero-indexing
#	shiftZero  <- shiftZero[0:(length(shiftZero)-4)]
#	shiftOne   <- shiftOne[0:(length(shiftOne)-3)]
#	shiftTwo   <- shiftTwo[0:(length(shiftTwo)-2)]
#	shiftThree <- shiftThree[0:(length(shiftThree)-1)]
#	
#	initialSum <- shiftZero + shiftOne + shiftTwo + shiftThree + shiftFour
#	
#	aggFinal   <- initialSum[c(TRUE, FALSE, FALSE, FALSE, FALSE)]
#	
#	return(aggFinal)
#}

# TR deprecated 20 July, 2018, use groupAges() instead

##' aggregates split 0 & 1-4 age groups into a single 5 year age group
##' 
##' @description Creates a five year age group from split 0 & 1-4 year age groups. 
##' @details Sums 0 & 1-4 age groups and outputs new 5 year age group vector.
##' 
##' @param Value numeric vector of population age groups that includes 0, 1-4, and 5 year ages.
##' 
##' @export 
##' @examples 
##' MalePop <- seq(1,100)
##' convertSplitTo5Year(MalePop)
#
#convertSplitTo5Year <- function(Value){
#	output <- rep(0, length(Value))
#	
#	intermediate1 <- Value[1]
#	intermediate2 <- Value[2]
#	
#	intermediateValue <- Value[-1]
#	intermediateValue[1] <- intermediate1 + intermediate2
#	
#	output <- data.frame(intermediateValue)
#	row.names(output) <- seq(0, 5*length(output[,1])-1, by = 5)
#	
#	return(output)
#}
