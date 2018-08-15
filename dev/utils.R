
# Author: tim
###############################################################################

#' shift a vector left or right
#' 
#' @description a simple function to move the elements of a vector to earlier or 
#' later positions depending on the value of \code{shift}. Auxiliary to other functions,
#'  rather multipurpose.

#' @param x a vector
#' @param shift an integer value \code{< length(x)}
#' @param fill what value do we put in the newly created positions? Possibly values like \code{FALSE}, \code{NA}, or  \code{0}.

#' @details Nothing fancy here. This is used for example in \code{Noumbissi()} to match denominator ranges to numerator positions using logical vectors.
#' @return the vector x, shifted left or right.
#' @export
shift.vector <- function(x,shift = 0, fill = FALSE){
	n <- length(x)
	
	if (shift > 0){
		x <- c(rep(fill,shift),x[-((n-shift+1):n)])
	}
	if (shift < 0){
		x <- c(x[-(1:abs(shift))],rep(fill,abs(shift)))
	}
	x
}



#' a simple centered moving average function

#' @description This function is defined based on a code chunk found here \url{https://stackoverflow.com/questions/743812/calculating-moving-average}. This is a centered moving average of arbitrary width. 
#' 
#' @param x numeric vector to produce a moving average of
#' @param n integer the width of the moving avarage. Default 5 time steps (years)
#' 
#' @return numeric vector of same length as \code{x}.
#' 
#' @details \code{NA} values
#' are used as padding on the left and right to make the returned vector of equal length.
#' 
#' @export
#' @examples 
#' x  <- runif(100)
#' xx <- cumsum(x)
#' \dontrun{
#' plot(x)
#' lines(ma(x))
#' lines(ma(x),9)
#' }


ma <- function(x,n=5){
	as.vector(stats::filter(x, rep(1 / n, n), sides = 2))
}


#' rescale a vector proportionally to a new sum

#' @description THis is a frequently needed operation, so it's being added as a general utility.
#' 
#' @param x numeric vector
#' @param scale what you want the vector to sum to. Default 1.
#' 
#' @details For a distribution, use \code{scale = 1}. For percentages, use \code{scale = 100}, etc.
#' 
#' @return the vector, rescaled
#' @examples 
#' x <- runif(10)
#' sum(x)
#' xx <- rescale.vector(x,100)
#' sum(xx)
#' @export

rescale.vector <- function(x, scale = 1){
	scale * x / sum(x, na.rm = TRUE)
}

#' @title determine whether a year is a leap year. 
#' 
#' @description In order to remove \code{lubridate} dependency, we self-detect leap years and adjust February accordingly. Code inherited from HMD.
#' 
#' @param Year integer of year to query
#' 
#' @return logical is the Year a leap year or not
#' 
#' @export
#' @author Carl Boe

is_LeapYear <- function (Year){      # CB: mostly good algorithm from wikipedia
	ifelse(
			( (Year %% 4) == 0  &  (Year %% 100) != 0   ) | ( (Year %% 400) == 0 ),
			TRUE, FALSE )
}

#' @title determine the proportion of a year passed as of a particular date
#' 
#' @description The fraction returned by this is used e.g. for intercensal estimates.
#' 
#' @param Year 4-digit year (string or integer)
#' @param Month month digits (string or integer, 1 or 2 characters)
#' @param Day Day of month digits (string or integer, 1 or 2 characters)
#' @param detect.mid.year logical. if \code{TRUE}, June 30 or July 1 will always return .5.
#' @param detect.start.end logical. default \code{TRUE}. Should Jan 1 always be 0 and Dec 31 always be 1?
#' @details Code inherited from HMD, slightly modified to remove matlab legacy bit.
#' @export
#' @examples
#' ypart(2001,2,14) # general use
#' ypart(2001,6,30) # mid year options, default detection
#' ypart(2001,7,1)  # also
#' ypart(2000,6,30) # no change for leap year, still detected as mid year
#' ypart(2000,6,30, FALSE) # exact measure
#' ypart(2000,7,1, FALSE)  # July 1 leap year
#' ypart(2001,7,1, FALSE)  # July 1 not leap year
#' ypart(2002,12,31, detect.start.end = FALSE) # assumes end of day by default.
#' ypart(2002,1,1, detect.start.end = FALSE) # end of day year fraction
#' ypart(2002,1,1, detect.start.end = TRUE)  # assume very begining of year


ypart <- function(Year, Month, Day, detect.mid.year = TRUE, detect.start.end = TRUE){
	M <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)

	if (detect.mid.year){
		.d <- as.integer(Day)
		.m <- as.integer(Month)
		if ((.d == 30 & .m == 6) | (.d == 1 & .m == 7)){
			return(.5)
		}
	}
	
	if (detect.start.end){
		.d <- as.integer(Day)
		.m <- as.integer(Month)
		if (.d == 1 & .m == 1){
			return(0)
		}
		if(.d == 31 & .m == 12){
			return(1)
		}
	}
	
	monthdur    <- diff(c(M,365))
	monthdur[2] <- monthdur[2] + is_LeapYear(Year)
	M           <- cumsum(monthdur) - 31
	return((M[Month] + Day) / sum(monthdur))
	
}


#' a convert date to decimal year fraction
#'  
#' @description Convert a character or date class to decimal, taking into account leap years. 

#' @details This makes use of two HMD functions, \code{ypart()}, and \code{is_LeapYear()} to compute. If the date is numeric, it is returned as such. If it is \code{"character"}, we try to coerce to \code{"Date"} class, ergo, it is best to specify a character string in an unambiguous \code{"YYYY-MM-DD"} format. If \code{date} is given in a \code{"Date"} class it is dealt with accordingly.
#' @param date either a \code{Date} class object or an unambiguous character string in the format \code{"YYYY-MM-DD"}.
#' 
#' @return numeric expression of the date, year plus the fraction of the year passed as of the date.
#' @export
dec.date <- function(date){
	if (class(date) == "numeric"){
		return(date)
	}
	if (class(date) == "character"){
		date <- as.Date(date)
	}
	day 	<- as.numeric(format(date,'%d'))
	month 	<- as.numeric(format(date,'%m'))
	year 	<- as.numeric(format(date,'%Y'))
	frac    <- ypart(Year = year, 
		                  Month = month, 
		                  Day = day,
		                  detect.mid.year = TRUE,
						  detect.start.end = TRUE)
	year + frac
}


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


#' take consecutive ratios of a vector
#' 
#' @description This can be used, for example to take survival ratios. 
#' @details behavior similar to \code{diff()}, in that returned vector is \code{k} elements
#' shorter than the given vector \code{fx}.
#' 
#' @param fx numeric vector of \code{length > k}
#' @param k integer the size of the lag in elements of \code{fx}.
#' 
#' @export 
#' @examples 
#' fx <- 1:10
#' ratx(fx)
#' ratx(fx,-1)
#' ratx(fx,0)
ratx <- function(fx, k = 1){
	k    <- as.integer(k)
	N    <- length(fx)
	m    <- N - abs(k)
	if (k > 0){
		fx <- fx[-c(1:k)] / fx[1:m]
	}
	if (k < 0){
		fx <- fx[1:m] / fx[-c(1:abs(k))] 
	}
	fx
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

#' infer abrdiged age group widths based on the length of a vector
#' 
#' @description This function is an auxiliary used strictly by top level functions where it is 
#' guaranteed that age groups are standard abridged age groups.
#' 
#' @param vec any vector, presumably a count, rate, or similar.
#' @param OAG logical (default \code{FALSE}). Is the final age group open?
#' 
#' @details This will give erroneous results if ages are anything other than standard abridged ages groups. 
#' If the final age group is open, the interval with is defined as \code{NA}. \code{Inf} or \code{-1} would have 
#' also been a good choice, but we went with \code{NA}.
#' 
#'  
#' @return an integer vector of \code{length(vec)} indicating the width of the abridged age group that each 
#' vector element corresponds to.
#' 
#' @export
#' @examples
#' vec <- runif(20)
#' inferAgeAbr(vec)
#' inferAgeAbr(vec, OAG = TRUE)
inferAgeAbr <- function(vec, OAG = FALSE){
	ageint    <- rep(5,length(vec))
	ageint[1] <- 1
	ageint[2] <- 4
	if (OAG){
		ageint[length(ageint)] <- NA
	}
	ageint
}


#' aggregates single year age groups into 5 year age groups
#' 
#' @description Creates five year age groups from single year ages. 
#' @details Sums five year age intervals
#' 
#' @param Value numeric vector of single year age groups.
#' 
#' @export 
#' @examples 
#' MalePop <- seq(1,100)
#' convertSingleTo5Year(MalePop)

convertSingleTo5Year <- function(Value){
  shiftZero  <- Value
  shiftOne   <- Value[-1]
  shiftTwo   <- shiftOne[-1]
  shiftThree <- shiftTwo[-1]
  shiftFour  <- shiftThree[-1]
  
  # TR: not sure what to make of zero-indexing
  shiftZero  <- shiftZero[0:(length(shiftZero)-4)]
  shiftOne   <- shiftOne[0:(length(shiftOne)-3)]
  shiftTwo   <- shiftTwo[0:(length(shiftTwo)-2)]
  shiftThree <- shiftThree[0:(length(shiftThree)-1)]
  
  initialSum <- shiftZero + shiftOne + shiftTwo + shiftThree + shiftFour
  
  aggFinal   <- initialSum[c(TRUE, FALSE, FALSE, FALSE, FALSE)]
  
  return(aggFinal)
}

#' aggregates split 0 & 1-4 age groups into a single 5 year age group
#' 
#' @description Creates a five year age group from split 0 & 1-4 year age groups. 
#' @details Sums 0 & 1-4 age groups and outputs new 5 year age group vector.
#' 
#' @param Value numeric vector of population age groups that includes 0, 1-4, and 5 year ages.
#' 
#' @export 
#' @examples 
#' MalePop <- seq(1,100)
#' convertSplitTo5Year(MalePop)

convertSplitTo5Year <- function(Value){
  output <- rep(0, length(Value))
  
  intermediate1 <- Value[1]
  intermediate2 <- Value[2]
  
  intermediateValue <- Value[-1]
  intermediateValue[1] <- intermediate1 + intermediate2
  
  output <- data.frame(intermediateValue)
  row.names(output) <- seq(0, 5*length(output[,1])-1, by = 5)
  
  return(output)
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
#' @return vector of counts in N-year age groups
#' 
#' @details If you shift the groupings, then the first age groups may have a negative lower bound 
#' (for example of -5). These counts would be discarded for the osculatory version of Sprague smoothing,
#' for example, but they are preserved in this function. The important thing to know is that if you shift 
#'  the groups, the first and last groups won't be N years wide. For example if \code{shiftdown} is 1, the first age group is 4-ages wide. The ages themselves are not returned, 
#' but they are the name attribute of the output count vector. Note this will also correctly group abridged ages
#' into equal 5-year age groups if the \code{Age} argument is explicitly given.
#' @export
#' @examples
#' India1991males <- c(9544406,7471790,11590109,11881844,11872503,12968350,11993151,10033918
#' 		,14312222,8111523,15311047,6861510,13305117,7454575,9015381,10325432,
#' 9055588,5519173,12546779,4784102,13365429,4630254,9595545,4727963
#' 		,5195032,15061479,5467392,4011539,8033850,1972327,17396266,1647397,
#' 6539557,2233521,2101024,16768198,3211834,1923169,4472854,1182245
#' 		,15874081,1017752,3673865,1247304,1029243,12619050,1499847,1250321,
#' 2862148,723195,12396632,733501,2186678,777379,810700,7298270
#' 		,1116032,650402,1465209,411834,9478824,429296,1190060,446290,362767,
#' 4998209,388753,334629,593906,178133,4560342,179460
#' 		,481230,159087,155831,1606147,166763,93569,182238,53567,1715697,
#' 127486,150782,52332,48664,456387,46978,34448
#' 		,44015,19172,329149,48004,28574,9200,7003,75195,13140,5889,18915,21221,72373)
#' Age <- 0:100
#' groupAges(India1991males, N = 5)
#' groupAges(India1991males, N = 5, shiftdown = 1)
#' groupAges(India1991males, N = 5, shiftdown = 2)
#' groupAges(India1991males, N = 5, shiftdown = 3)
#' groupAges(India1991males, N = 5, shiftdown = 4)
groupAges <- function(Value, 
		              Age = 1:length(Value) - 1, 
					  N = 5, 
					  shiftdown = 0, 
					  AgeN){
	if (missing(AgeN)){
		AgeN <- calcAgeN(Age = Age, N = N, shiftdown = shiftdown)
	}
	tapply(Value, AgeN, sum)
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
#' test.single(Age)  # TRUE
#' test.single(Age2) # FALSE repeated, can't tell.
#' test.single(Age3) # FALSE not single ages
#' test.single(Age4) # FALSE not single ages
test.single <- function(Age){
	all(diff(sort(Age)) == 1)
}

#' Convert arbitrary age groupings into single years of age
#' 
#' @description Splits aggregate counts for a vector of age groups of a common width into single year age groups.
#' 
#' @param Value numeric vector of age group counts
#' @param Age numeric vector of ages corresponding to the lower integer bound of the age range.
#' @param OAG boolean argument that determines whether the final age group (assumed open ended) is kept as it is or has the same length as the rest of the age groups. Default is FALSE, i.e. use the same length for the final age group.
#' 
#' @return numeric vector of counts for single year age groups.
#' 
#' @details Assumes that all age groups are equal width. (Default is 5-year age groups.) Also, assumes that the population is uniformly distributed across each age interval. If there is a split under 5 age group (0 and 1-4 age groups), use function convertSplitTo5Year to consolidate before using this function. The default setting is to assume that the final age group is the same width as the other age groups.
#' 
#' @export
#' @examples 
#' MalePop <- c(9544406,7471790,11590109,11881844,11872503,12968350,
#'              11993151,10033918,14312222,8111523,15311047,6861510,13305117,7454575,
#'              9015381,10325432,9055588,5519173)
#' Ages <- seq(0,85, by=5)
#' splitToSingleAges(MalePop, Ages)
#' splitToSingleAges(MalePop, Ages, OAG = TRUE)
splitToSingleAges <- function(Value, Age, OAG = FALSE){
  ageMax <- max(Age)             # the lower bound of the largest age group
  N      <- length(Value)        # number of age groups from the census.
  M      <- ageMax / (N - 1)     # length of each age group from the census.
  
  ageGroupBirths   <- Value / M  # vector of the number of births in a single year for each age group assuming uniformity.
  
  if (OAG){
      ageGroupBirths <- ageGroupBirths[0:(length(ageGroupBirths)-1)] # remove final age group
      singleAgeGroupBirths <- rep(ageGroupBirths, each = M) # vector of the single year births for all ages except final age group
      singleAgeGroupBirths[length(singleAgeGroupBirths)+1] <- Value[length(Value)]
  }
  else{
      singleAgeGroupBirths  <- rep(ageGroupBirths, each = M)  # vector of the single year births for all ages
  }
  
  return(singleAgeGroupBirths)
}

#' Wrapper to provide a single location to reference all model life tables.
#' @description Still in the works.
#' @param ModelName string naming the life table to return. Can be "coale-demeny west".
#' @param Sex string indicating which sex should be returned. Can be either "M" or "F".
#' @return list of life tables
#' @details More model families can easily be added.
#' @importFrom demogR cdmltw
#' @export

getModelLifeTable <- function(ModelName, Sex){
	Sex <- toupper(substr(Sex, 1, 1))
	
	stopifnot(Sex %in% c("M", "F"))
	stopifnot(ModelName == "coale-demeny west")
 
    outputLT <- demogR::cdmltw(sex = Sex)
  
  
  return(outputLT)
}



# TODO in progress

#' rescale counts in age groups to match counts in different age groups
#' @description This method rescales a vector of counts in arbitrary (integer) age groups to approximate a vector of counts in a potentially different age grouping. Common use cases will be to scale single ages (whose age pattern we wish to roughly maintain) to sum to abridged or 5-year age groups from another source. The counts to be rescaled could potentially be in any grouping (see example). 
#' @details If the final age group is open, define its age interval as 1.
#' 
#' Presently the intermediate splitting function assumes that counts inside the age groups of population 1 are uniformly distributed, although this may be relaxed if other methods become available whose behavior matches that of \code{splitUniform()}. \code{splitMono()} will be modified soon to be applicable here. 
#' 
#' The method works by first spl
#' @examples
#' # just to make a point about arbitrary integer age widths in both pop1 and pop2
#' # note if pop1 is in single ages and pop2 is in groups things work much cleaner.
#' set.seed(3)
#' AgeIntRandom <- sample(1:5, size = 15,replace = TRUE)
#' AgeInt5      <- rep(5, 9)
#' original     <- runif(45, min = 0, max = 100)
#' pop1         <- groupAges(original, 0:45, AgeN = int2ageN(AgeIntRandom, FALSE))
#' pop2         <- groupAges(original, 0:45, AgeN = int2ageN(AgeInt5, FALSE))
#' # inflate (in this case) pop2
#' perturb      <- runif(length(pop2), min = 1.05, max = 1.2)
#' 
#' pop2         <- pop2 * perturb
#' 
#' # a recursively constrained solution
#' (pop1resc <- rescaleAgeGroups(Value1 = pop1, 
#' 					AgeInt1 = AgeIntRandom,
#' 					Value2 = pop2, 
#' 					AgeInt2 = AgeInt5, 
#' 					splitfun = splitUniform, 
#' 					res = 1e-3,
#'                  recursive = TRUE))
#' # a single pass adjustment (no recursion)
#' (pop1resc1 <- rescaleAgeGroups(Value1 = pop1, 
#' 					AgeInt1 = AgeIntRandom,
#' 					Value2 = pop2, 
#' 					AgeInt2 = AgeInt5, 
#' 					splitfun = splitUniform, 
#' 					res = 1e-3,
#' 					recursive = FALSE))
#' pop1resc / pop1
#' perturb
#' 
#' 
#' \dontshow{
#' 	newN        <- splitfun(pop1resc, AgeInt = AgeInt1)
#' 	check       <- groupAges(newN, AgeS, AgeN = AgeN2)
#' 	stopifnot(all(abs(check - pop2) < 1e-3))
#' # and the non-recursive one:
#' 	newN        <- splitfun(pop1resc1, AgeInt = AgeInt1)
#' 	check       <- groupAges(newN, AgeS, AgeN = AgeN2)
#' 	stopifnot(all(abs(check - pop2) < 1e-3))
#' }
#' 
#' \dontrun{
#' 	
#' # show before / after
#' 	plot(NULL,xlim=c(0,45),ylim=c(0,2),main = "Different (but integer) intervals",
#' 			xlab = "Age", ylab = "", axes = FALSE)
#' 	x1 <- c(0,cumsum(AgeIntRandom))
#' 	rect(x1[-length(x1)],1,x1[-1],2,col = gray(.8), border = "white")
#' 	x2 <- c(0,cumsum(AgeInt5))
#' 	rect(x2[-length(x2)],0,x2[-1],1,col = "palegreen1", border = "white")
#' 	text(23,1.5,"Original (arbitrary grouping)",font = 2, cex=1.5)
#' 	text(23,.5,"Standard to rescale to (arbitrary grouping)",font = 2, cex=1.5)
#' 	axis(1)
#' 	
#' # adjustment factors:
#' 	plot(int2age(AgeInt5), perturb, ylim = c(0, 2))
#' 	points(int2age(AgeIntRandom), pop1resc / pop1, pch = 16)
#' # non-recursive is less disruptive for uniform
#' 	points(int2age(AgeIntRandom), pop1resc1 / pop1, pch = 16, col = "blue")
#' 	
#' # show before / after under uniform (in pop1) assumption.
#' 	plot(NULL, xlim=c(0,45),ylim = c(0,150), main = "Uniform constraint")
#' 	lines(0:44, splitUniform(pop1, AgeInt = AgeIntRandom, OAG = FALSE), col = "red")
#' 	lines(0:44, splitUniform(pop2, AgeInt = AgeInt5, OAG = FALSE), col = "blue")
#' 	lines(0:44, splitUniform(pop1resc, AgeInt = AgeIntRandom, OAG = FALSE), col = "orange", lty = 2, lwd = 2)
#' 	lines(0:44, splitUniform(pop1resc1, AgeInt = AgeIntRandom, OAG = FALSE), col = "magenta", lty = 2, lwd = 2)
#' 	legend("topright",lty=c(1,1,2,2),col=c("red","blue","orange","magenta"),lwd=c(1,1,2,2),
#' 			legend = c("Original N1","Prior N2","Rescaled N1 recursive","Rescaled N1 1 pass"))
#' 	
#' }
#' 

rescaleAgeGroups <- function(
		Value1, 
		AgeInt1, 
		Value2, 
		AgeInt2, 
		splitfun = splitUniform, 
		recursive = FALSE,
		res = 1e-3){
	N1          <- length(Value1)
	# ages must cover same span
	stopifnot(sum(AgeInt1) == sum(AgeInt2))
	
	Age1        <- int2age(AgeInt1)
	Age2        <- int2age(AgeInt2)
	
	stopifnot(N1 == length(Age1))

	AgeN        <- rep(Age2, times = AgeInt2)
	
	# scale from single.
	ValueS      <- splitfun(Value1, AgeInt = AgeInt1)
	# right now splitMono() doesn't have AgeInt, so does not create the right spread.
	# comparison forthcoming.
	AgeS        <- names2age(ValueS)
	
	#splitMono(Value = Value1, Age1, F)
	#plot(splitMono(Value = Value1, Age1, T))

	#
	AgeN2       <- rep(Age2, times = AgeInt2)
	beforeN     <- groupAges(ValueS, AgeS, AgeN = AgeN2)
	
	beforeNint  <- rep(beforeN, times = AgeInt2)
	afterNint   <- rep(Value2, times = AgeInt2)
	ratio       <- afterNint / beforeNint
	SRescale    <- ValueS * ratio
	
	
	
	# group back to original, in case these weren't single
	AgeN1       <- rep(Age1, times = AgeInt1)
	out         <- groupAges(SRescale, AgeS, AgeN = AgeN1)
	
	# check for recursion
	newN        <- splitfun(out, AgeInt = AgeInt1)
	check       <- groupAges(newN, AgeS, AgeN = AgeN2)
	if (max(abs(check-pop2)) < res | !recursive){
		return(out)
	} else {
		rescaleAgeGroups(
				Value1 = out,  
				AgeInt1 = AgeInt1, 
				Value2 = Value2, 
				AgeInt2 = AgeInt2,
				splitfun = splitfun,
				res = res,
				recursive = recursive)
	}
}



