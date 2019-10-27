
# Author: tim
###############################################################################
# contains utilities having to do with age groups.
# - simple operations like detecting, grouping, and so forth.
# not complex things like ungrouping ages.

#' Trim two age vectors to matching N-year age groups.
#' 
#' @description Determining which N-year (e.g. 5-year) age groups two vectors have in common
#' is helpful for matching vector lengths, and questions of alignment. Used as a utility throughout.

#' @param Age1 integer. Vector of first age groups (lower bounds).
#' @param Age2 integer. Vector of second age groups (lower bounds).
#' @param N integer. Target age group interval (just one number).
#' @param consecutive logical.  Whether or not to throw error if resulting age groups not consecutive. Default \code{TRUE}.
#' @param ageMin integer. Optional lower age bound for output.
#' @param ageMax integer. Optional upper age bound for output.
#' 
#' @return Integer vector of the N-year age groups present in both \code{Age1} and \code{Age2}.
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

#' Calculate which large age group single ages belong to.
#' 
#' @description Assign single ages to age groups of equal and arbitrary width, and also optionally shifted.
#' 
#' @param Age integer. Vector of single ages (lower bound).
#' @param N integer. Desired width of resulting age groups.
#' @param shiftdown integer. Move the grouping down by one or more single ages. Optional argument.
#' 
#' @details If you shift the groupings, then the first age groups may have a negative lower bound 
#' (for example of -5). These counts would be discarded for the osculatory version of Sprague smoothing,
#' for example, but they are preserved in this function. The important thing to know is that if you shift 
#'  the groups, the first and last groups won't be N years wide. For example if \code{shiftdown} is 1, 
#' the first age group is 4-ages wide. 
#'  
#' @return An integer vector of \code{length(Age)} indicating the age group that each single age belongs to.
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

#' repeat age lower bounds once for each single age
#' @description This is a frequent grouping situation. For a given vector of lower age bounds, repeat each value N times, where N is the width of the corresponding age interval. Age intervals are in this case calculated from the original Age vector.
#' @details If {OAG = TRUE} then the last value is not repeated, otherwise the final age interval is assumed to be the same width as the penultimate. Here intervals do not need to be of uniform width.
#' @param Age integer. Vector of lower age bounds.
#' @param OAG logical. Whether or not the final age group open. Default \code{FALSE}. See details
#' @return integer vector of elements of \code{Age} repeated once for each implied single age.
#' @export
#' @examples 
#' age1 <- seq(0,100,by=5)
#' (ageN1 <- age2ageN(age1, OAG = FALSE))
#' (ageN2 <- age2ageN(age1, OAG = TRUE))
age2ageN <- function(Age, OAG = FALSE){
	rep(Age, times = age2int(Age, OAG = OAG, OAvalue = 1))
}

#' repeat age lower bounds once for each single age
#' @description This is a frequent grouping situation. For a given vector of lower age bounds, as implied by \code{AgeInt}, repeat each value N times, where N is the width of the corresponding age interval. Age intervals are in this case given, and age lower bound values are calculated from \code{AgeInt} vector.
#' @details If {OAG = TRUE} then the last value is given just once, irrespective of the final value of \code{AgeInt}, otherwise the final age interval is repeated \code{AgeInt[length(AgeInt)]} times. Here intervals do not need to be of uniform width.
#' @param AgeInt integer or numeric. Vector of age intervals.
#' @param OAG logical. Whether or not the final age group open. Default \code{FALSE}. See details
#' @return integer vector of elements of \code{Age} repeated once for each implied single age.
#' @export
#' @examples 
#' int5 <- rep(5,21)
#' (ageN1 <- int2ageN(int5, OAG = FALSE))
#' (ageN2 <- int2ageN(int5, OAG = TRUE))
int2ageN <- function(AgeInt, OAG){
	if (OAG){
		AgeInt[length(AgeInt)] <- 1
	}
	rep(int2age(AgeInt), times = AgeInt)
}

#' Calculate which abridged age group single ages belong to.
#' 
#' @description Assign single ages to 5-year abridged age groups. That means that age 0 is kept as a single age,
#' ages 1-4 are grouped together as abridged age 1, and thereafter 5-year age groups are used.
#' 
#' @param Age integer. Vector of single ages (lower bound).
#' 
#' @details In the case that the single \code{Age} vector starts at some age higher than 4, 
#' this is just the same as \code{calcAgeN(Age,5,0)}.
#'  
#' @return An integer vector of \code{length(Age)} indicating the abridged age group that each single age belongs to.
#' 
#' @export
#' @examples
#' Age <- 0:70
#' calcAgeAbr(Age)
#' calcAgeN(Age,5,0)
calcAgeAbr <- function(Age){
	stopifnot(is.integer(Age))
	Abr               <- Age - Age %% 5
	Abr[Age %in% 1:4] <- 1
	Abr
}

#' Infer abridged age groups widths.
#' 
#' @description This function is an auxiliary used by top level functions where it is 
#' guaranteed that age groups are standard abridged age groups. If \code{Age} is specified,
#' this will work as well for other age groupings.
#' 
#' @param Age integer. Vector of lower bound of each age group.
#' @param vec Any vector, presumably a count, rate, or similar.
#' @param OAG logical. Whether or not the final age group open. Default \code{FALSE}. 
#' @param OAvalue numeric or integer. The value to use for the final age interval if \code{OAG = TRUE}. Default \code{NA}.
#' 
#' @details If based solely on the length of a vector, this will give erroneous results if ages 
#' are anything other than standard abridged ages groups. If the final age group is open, the 
#' interval width is defined as \code{NA}. \code{Inf} or \code{-1} would have 
#' also been a good choice, but we went with \code{NA}.
#' 
#'  
#' @return An integer vector of \code{length(vec)} indicating the width of the abridged age group that each 
#' vector element corresponds to.
#' 
#' @export
#' @examples
#' vec <- runif(20)
#' inferAgeIntAbr(vec = vec)
#' inferAgeIntAbr(vec = vec, OAG = TRUE)
inferAgeIntAbr <- function(Age, vec, OAG = FALSE, OAvalue = NA){
	
	# Age is preferred (lower bounds)
	if (!missing(Age)){
		ageint <- age2int(Age = Age, 
				          OAG = OAG, 
						  OAvalue = OAvalue)
	}
	# otherwise length of vector will do
	if (missing(Age) & !missing(vec)){
		stopifnot(length(vec)>5)
		ageint    <- rep(5,length(vec))
		ageint[1] <- 1
		ageint[2] <- 4
		if (OAG){
			ageint[length(ageint)] <- OAvalue
		}
	}
	
	ageint
}

#' Determine abridged ages up to a given maximum age group.
#' @description Produce standard abridged age groups (lower bounds) up to a specified maximum age group. 
#' @details If the highest age group is not evenly divisible by 5 then age classes only go up to its 5-year lower bound.
#' @param ageMax integer. Default 80.
#' @return integer. Vector of ages, e.g. \code{c(0,1,5,10,15,...)}.
#' @export 
#' @examples 
#' maxA2abridged(80)
#' all(maxA2abridged(100) == maxA2abridged(102))
maxA2abridged <- function(ageMax = 80){
	sort(unique(calcAgeAbr(0:ageMax)))
}



#' Infer lower age bounds from age class intervals.
#' @description Determine lower bounds of age classes based on a vector of age intervals and a starting age.
#' 
#' @param AgeInt integer or numeric. Vector of age intervals.
#' @param ageMin integer. The lowest age, default 0.
#' @export 
#' @return Age vector of same length as \code{AgeInt}.
#' @examples 
#' AgeInt <- c(1,4,rep(5,17))
#' int2age(AgeInt)
int2age <- function(AgeInt, ageMin = 0){
  n <- length(AgeInt)
  # if final AgeInt is NA, then assume it's OAG,
  # count as zero for this calc
  if (is.na(AgeInt[n])){
    AgeInt[n] <- 0
  }
	cumsum(AgeInt) - AgeInt + ageMin
}

#' Infer age class intervals from lower age bounds.
#' @description Determine age class intervals based on a vector of age class lower bounds.
#' @details If the final age group is open, it is given a value of \code{NA} by default, or else a user-determined value. 
#' If the final age group is closed, it is assumed to be equal to the next-lower interval. If the final age interval is 
#' known and not equal to the next lowest interval, specify \code{OAG = TRUE} and assign its value to \code{OAvalue}.
#' @param Age integer or numeric. Vector of lower age group bounds . 
#' @param OAG logical. Whether or not the final age group is open. Default \code{TRUE}.
#' @param OAvalue numeric or integer. The value to use for the final age interval if \code{OAG = TRUE}. Default \code{NA}.
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


#' Group single ages into equal age groups of arbitrary width
#' 
#' @description This can be useful to check constrained sums, or as an intermediate step for smoothing.
#' 
#' @param Value numeric. Vector of single age counts.
#' @param Age integer. Vector of lower bounds of single age groups.
#' @param N integer. The desired width of resulting age groups. Default 5.
#' @param shiftdown integer. Optionally shift age groupings down by single ages. Default 0.
#' @param AgeN integer vector, otherwise calculated using \code{calcAgeN()}. Optional argument.
#' @param OAnew integer. Value of lower bound of new open age group.
#' @return Vector of counts in N-year age groups.
#' 
#' @details If you shift the groupings, then the first age groups may have a negative lower bound 
#' (for example of -5). These counts would be discarded for the osculatory version of Sprague smoothing,
#' for example, but they are preserved in this function. The important thing to know is that if you shift 
#'  the groups, the first and last groups will not be N years wide. For example if \code{shiftdown} is 1, the first age group is 4-ages wide. The ages themselves are not returned, 
#' but they are the name attribute of the output count vector. Note this will also correctly group abridged ages
#' into equal 5-year age groups if the \code{Age} argument is explicitly given. \code{OAnew} (optional) must be less than or equal to \code{max(Age)} to have any effect.
#' @export
#' @examples
#'  Age <- 0:100
#'  groupAges(pop1m_ind, N = 5)
#'  groupAges(pop1m_ind, N = 5, shiftdown = 1)
#'  groupAges(pop1m_ind, N = 5, shiftdown = 2)
#'  groupAges(pop1m_ind, N = 5, shiftdown = 3)
#'  groupAges(pop1m_ind, N = 5, shiftdown = 4)
#'  groupAges(pop1m_ind, N = 5, OAnew = 80) 

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

#' Group down to a new open age class.
#' 
#' @description This simple utility lowers the open age group. It only returns the input value vector, not the age vector.
#' @param Value numeric. Vector of counts.
#' @param Age integer. Vector of age classes.
#' @param OAnew The desired open age group.
#' @export 
#' @return Value vector potentially of reduced length up to OAG.

groupOAG <- function(Value, Age, OAnew){
	stopifnot(OAnew %in% Age)
	N        <- length(Value[Age <= OAnew])
	Value[N] <- sum(Value[Age >= OAnew])
	Value    <- Value[1:N]
	Value
}


#' check for coherence within Age and between Age and AgeInt
#' 
#' @description A few checks are carried out to test if \code{Age} is internally consistent, that \code{OAG} is consistent with \code{AgeInt}, and that \code{Age} and \code{AgeInt} are consistent with one another. For \code{Age} to be internally consistent, we cannot have redundant values, and values must be sequential. 
#' @details If \code{OAG} is \code{TRUE} then \code{AgeInt} must be coded as \code{NA}. If \code{Age} is not sorted then we sort both \code{Age} and \code{AgeInt}, assuming that they are in matched order. This isn't incoherence per se, but a message is returned to the console.
#' 
#' @param Age integer vector of single ages (lower bound)
#' @param AgeInt integer vector. Age interval widths
#' @param OAG logical (default \code{TRUE}). Is the final age group open?
#' @return logical. \code{TRUE} if the arguments are considered consistent.
#' @export
#' @examples 
#' Age    <- 0:99
#' AgeInt <- rep(1, 100)
#' # closed, sequential, non-redundant ages, any easy yes:
#' is_age_coherent(Age = Age, AgeInt = AgeInt, OAG = FALSE)     # TRUE
#' 
#' # incorrectly coded OAG
#' is_age_coherent(Age = Age, AgeInt = AgeInt, OAG = TRUE)      # FALSE   
#' 
#' # correctly coded OAG
#' AgeInt[100] <- NA
#' is_age_coherent(Age = Age, AgeInt = AgeInt, OAG = TRUE)      # TRUE
#' 
#' # correct unordered, but this isn't incoherence per se.
#' # watch out though!
#' aaoo <- order(sample(Age, 100, replace = FALSE))
#' is_age_coherent(Age[aaoo], AgeInt = AgeInt[aaoo], OAG = TRUE) # TRUE
#' 
#' # check redundancy
#' AgeRed    <- c(0:100,70)
#' AgeInt    <- c(rep(1, 100), NA, NA)
#' ao        <- order(AgeRed)
#' AgeRed    <- AgeRed[ao]
#' AgeIntRed <- AgeInt[ao]
#' is_age_coherent(AgeRed, AgeInt, TRUE)  # FALSE
is_age_coherent <- function(Age, AgeInt, OAG = TRUE){
  n       <- length(Age)
  stopifnot(length(AgeInt) == n)
  
  if (!is_age_sequential(Age)){
    message("Age isn't sorted. Sorting Age and AgeInt\nunder strong assumption that they are in matched order.\n")
    ao     <- order(Age)
    Age    <- Age[ao]
    AgeInt <- AgeInt[ao]
  }
  
  # check redundancy (case of more than one open age)
  if (is_age_redundant(Age, AgeInt)){
    message("At least one age repeated, ergo Age is incoherent")
    return(FALSE)
  }
  
  # is Age implied by AgeInt same?
  Age2    <- int2age(AgeInt, ageMin = min(Age))
  
  # better to compare 0s than NAs in next check
  if (OAG){
    if (!is.na(AgeInt[n])){
      message("The AgeInt value for OAG should be NA, but it is", AgeInt[n])
      return(FALSE)
    }
    AgeInt[n] <- 0
  }
  # is AgeInt implied by Age same?
  AgeInt2 <- age2int(Age, OAG = OAG, OAvalue = 0)
  
  # is everything coherent?
  out     <- all(Age == Age2) & all(AgeInt == AgeInt2)
  if (!out){
    message("Age and AgeInt don't appear to imply one another")
  }
  out
}



#' Logical checking of whether age classes appear single.
#' 
#' @description Check whether a vector of ages consists in single ages. This
#' makes sense sometimes when age intervals are not explicitly given.
#' 
#' @param Age integer. Vector of age classes.
#' 
#' @return logical. \code{TRUE} if detected as single ages, \code{FALSE} otherwise.
#' 
#' @details In cases where ages are indeed single, but repeated, this will still return \code{FALSE}. 
#' Therefore make sure that the age vector given refers to a single year of a single population.
#' @export
#' @examples
#' Age <- 0:99
#' Age2 <- c(0:10,0:10)
#' Age3 <- seq(0,80,by=5)
#' Age4 <- seq(0,10,by=.5)
#' is_single(Age)  # TRUE
#' is_single(Age2) # FALSE repeated, can't tell.
#' is_single(Age3) # FALSE not single ages
#' is_single(Age4) # FALSE not single ages
is_single <- function(Age){
	all(diff(sort(Age)) == 1)
}


#' is \code{Age} sorted sequentially?
#' 
#' @description Check if \code{Age} is sorted sequentially.
#' @details This does not check for redundancy.
#' @inheritParams is_age_coherent
#' @export
#' @return logical. Is the \code{Age} vector sorted sequentially?
is_age_sequential <- function(Age){
  all(Age == sort(Age))
}

#' check for redundant age specification
#' 
#' @description Ages are considered redundant if values for the underlying single ages are repeated. This might occur if there is an extra open age group below the final open age group. For example we have single ages 0 to 84, with an open age group of 85+, but the data also contain an open age group of 70+, leading to age 70 appearing twice. This will also detect 
#' @details Missing \code{AgeInt} are conservatively imputed with 20, which will most often trigger \code{FALSE} and thereby flag for further inspection. 
#' @inheritParams is_age_coherent
#' @return logical. Are there repeated values in the \code{Age} vector?
#' @export
#' @examples 
#' Age    <- c(0:100,70)
#' AgeInt <- c(rep(1,100),NA,NA)
#' # it doesn't matter if Age is sequential for this check
#' is_age_redundant(Age, AgeInt)
is_age_redundant <- function(Age, AgeInt){
  
  if (any(is.na(AgeInt))){
    AgeInt[is.na(AgeInt)] <- 20
  }
  AgeTo   <- Age + AgeInt - 1
 
  singles <- unlist(mapply(seq,Age,AgeTo))
  
  any(table(singles)>1)
}

#' Detect if a vector of lower age bounds is plausibly of abridged ages.
#' @description A logical utility that checks if a given vector is of the lower bounds of abridged age groups or not.
#' @param Age integer. Vector of lower age bounds.
#' @export
#' @return Logical \code{TRUE} if abridged, \code{FALSE} otherwise.

#' @examples
#' # as expected, TRUE
#' is_abridged(c(0,1,5,10,15,20,25))
#' # standard 5, not abridged, FALSE
#' is_abridged(c(0,5,10,15,20,25))
#' # plausible, TRUE
#' is_abridged(c(1,5,10,15,20,25))
#' # plausible, TRUE
#' is_abridged(c(5,10,15,20,25))
#' # 10 year age group not abridged, FALSE
#' is_abridged(c(0,1,5,10,15,25))
is_abridged <- function(Age){
	Age           <- as.integer(Age)
	ageMax        <- max(Age)
	ageMin        <- min(Age)
	abr_default   <- maxA2abridged(ageMax)
	abr_default   <- abr_default[abr_default >= ageMin]
	out           <- length(Age) == length(abr_default)
	if (out){
		out <- all(Age==abr_default)
	}
	out
}

#' detect ages from names of vector(s)
#' @description Often as a shorthand we pull lower age bounds from the names of a vector. This modulates that action, and allows for giving several vectors to check for names.
#' @details If more than one vector is given, names are pulled from the first available named vector. All vectors must be of the same length. If no names are available on any vector, then \code{NA}s are returned. This clearly won't work if the names on the vector are of something else.
#' @param ... one or more vectors of any class, which has/have a names attribute.
#' @return integer vector of ages (presumably).
#' @export 
#' @examples 
#' # create some vectors
#' nMx <- c(0.11621,0.02268,0.00409,0.00212,0.00295,0.00418,0.00509,0.00609,
#' 0.00714,0.00808,0.00971,0.0125,0.0175,0.02551,0.03809,0.05595,0.08098,
#' 0.15353,0.2557)
#' 
#' nqx <- c(0.10793,0.08554,0.02025,0.01053,0.01463,0.02071,0.02515,0.02999,
#' 0.03507, 0.03958,0.04742,0.0606,0.08381,0.11992,0.17391,0.2454,0.33672,
#' 0.54723,NA)
#' 
#' lx  <- c(100000,89207,81577,79924,79083,77925,76312,74393,72162,69631,66875,
#' 63704,59843,54828,48253,39861,30079,19951,9033)
#' age <- int2age(inferAgeIntAbr(vec=nMx))
#' names(nMx) <- age
#' names2age(nMx)
#' # or two vectors (only one has names here)
#' names2age(nMx,nqx)
#' # or more
#' names2age(nMx,nqx,lx)
#' # order doesn't matter
#' names2age(nqx,nMx, lx)
#' # multiple named, take first
#' names2age(nMx, nMx)
#' # NAs returned
#' names2age(lx)

names2age <- function(...){
	XL <- list(...)
	# if multiple vectors given must be same lengths
	if (length(XL) > 1){
		LL <- unlist(lapply(XL,function(x){
							length(x)
						}))
		stopifnot(all(diff(LL) == 0))
	}
	# which have names
	TF <- unlist(lapply(XL,function(x){
	      name_fun <- ifelse(is.null(dim(x)), names,function(y){dimnames(y)[[1]]})
				!is.null(name_fun(x))
			}))
	
	if (any(TF)){
		# takes names from first available
		x        <- XL[[which(TF)[1]]]
		name_fun <- ifelse(is.null(dim(x)), names,function(y){dimnames(y)[[1]]})
		Age      <- as.integer(name_fun(x))
	} else {
		x        <- XL[[1]]
		Age      <- rep(NA, length(x))
	}
	Age
}

#' rescale counts in age groups to match counts in different age groups
#' @description This method rescales a vector of counts in arbitrary (integer) age groups to approximate a vector of counts in a potentially different age grouping. Common use cases will be to scale single ages (whose age pattern we wish to roughly maintain) to sum to abridged or 5-year age groups from another source. The counts to be rescaled could potentially be in any grouping (see example). 
#' @details If the final age group is open, define its age interval as 1.
#' 
#' Presently the intermediate splitting function can either be \code{splitUniform()} or \code{splitMono()}. 
#' 
#' The method is an original contribution. It works by first splitting the counts of \code{Value1} to single ages using the assumptions of \code{splitfun()}. \code{Value1} is then rescaled such that were it re-grouped to match the age classes of \code{Value2} they would be identical. If \code{recursive = FALSE}, the single-age rescaled \code{Value1} data are returned regrouped to their original ages. If \code{recursive = TRUE}, the process is repeated until \code{Value1} is rescaled such that it could be split and regrouped to \code{Value2} using the same process a single time with no need for further rescaling. If age groups in \code{Value1} are very irregular, \code{recursive = TRUE} can induce noise (see example). If the age groups of \code{Value1} nest cleanly within the age groups of \code{Value2} then recursion is unnecessary. This is the case, for example, whenever \code{Value1} is in single ages and \code{Value2} is in grouped ages, which is likely the most common usage scenario.
#' @param Value1 numeric vector. A vector of demographic counts for population 1.
#' @param AgeInt1 integer vector. Age interval widths for population 1.
#' @param Value2 numeric vector. A vector of demographic counts for population 2.
#' @param AgeInt2 integer vector. Age interval widths for population 2.
#' @param splitfun function to use for splitting \code{pop1}. Presently on \code{splitUniform()} works.
#' @param recursive logical. Shall we repeat the split/regroup/rescale process until stable? See details. Default \code{FALSE}.
#' @param tol numeric. Default 1e-3. The numerical tolerance for the residual. Used to detect stability if \code{recursive = TRUE}.
#' @export
#' 
#' @examples
#' # just to make a point about arbitrary integer age widths in both pop1 and pop2
#' # note if pop1 is in single ages and pop2 is in groups things work much cleaner.
#' set.seed(3)
#' #set.seed(3)
#' #AgeIntRandom <- sample(1:5, size = 15,replace = TRUE)
#' AgeIntRandom <- c(1L, 5L, 2L, 2L, 4L, 4L, 1L, 2L, 3L, 4L, 3L, 3L, 3L, 3L, 5L)
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
#'                  recursive = TRUE))
#' # a single pass adjustment (no recursion)
#' (pop1resc1 <- rescaleAgeGroups(Value1 = pop1, 
#' 					AgeInt1 = AgeIntRandom,
#' 					Value2 = pop2, 
#' 					AgeInt2 = AgeInt5, 
#' 					splitfun = splitUniform, 
#' 					recursive = FALSE))
#' pop1resc / pop1
#' perturb
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
#' 	plot(NULL, xlim = c(0, 45), ylim = c(0, 150), main = "Uniform constraint")
#' 	lines(0:44, splitUniform(pop1, AgeInt = AgeIntRandom, OAG = FALSE), col = "red")
#' 	lines(0:44, splitUniform(pop2, AgeInt = AgeInt5, OAG = FALSE), col = "blue")
#' 	lines(0:44, splitUniform(pop1resc, AgeInt = AgeIntRandom, OAG = FALSE), 
#' col = "orange", lty = 2, lwd = 2)
#' 	lines(0:44, splitUniform(pop1resc1, AgeInt = AgeIntRandom, OAG = FALSE), 
#' col = "magenta", lty = 2, lwd = 2)
#' 	legend("topright",
#'         lty = c(1, 1, 2, 2),
#'         col = c("red", "blue", "orange", "magenta"),
#'         lwd = c(1, 1, 2, 2),
#' 		   legend = c("Original N1", "Prior N2",
#'                     "Rescaled N1 recursive", "Rescaled N1 1 pass"))
#' 	
#' }
#' 

rescaleAgeGroups <- function(
		Value1, 
		AgeInt1, 
		Value2, 
		AgeInt2, 
		splitfun = c(splitUniform,splitMono)[1], 
		recursive = FALSE,
		tol = 1e-3){
	N1          <- length(Value1)
	# ages must cover same span
	stopifnot(sum(AgeInt1) == sum(AgeInt2))
	
	Age1        <- int2age(AgeInt1)
	Age2        <- int2age(AgeInt2)
	
	stopifnot(N1 == length(Age1))
	
	AgeN        <- rep(Age2, times = AgeInt2)
	
	# step 1) split single. (innocuous if already single)
	ValueS      <- splitfun(Value1, AgeInt = AgeInt1)
	AgeS        <- names2age(ValueS)
	# right now splitMono() doesn't have AgeInt, so does not create the right spread.
	# comparison forthcoming.

	# step 2) regroup to groups of Value2 
	AgeN2       <- rep(Age2, times = AgeInt2)
	beforeN     <- groupAges(ValueS, AgeS, AgeN = AgeN2)
	
	# step 3) now repeat values of Value1 and Value2 for each single age
	# then rescale single age values.
	beforeNint  <- rep(beforeN, times = AgeInt2)
	afterNint   <- rep(Value2, times = AgeInt2)
	ratio       <- afterNint / beforeNint
	SRescale    <- ValueS * ratio
	
	# step 4) group back to original age classes
	AgeN1       <- rep(Age1, times = AgeInt1)
	out         <- groupAges(SRescale, AgeS, AgeN = AgeN1)
	
	# step 5a) if no recursion, return now
	if (!recursive){
		return(out)
	}
	# step 5b) otherwise continue, and only stop if
	# steps 1+2 from above would imply Value2 in a single pass.
	# Risky if an arbitrary splitting function is used...
	# equivalent of a while loop with no escape.
	newN        <- splitfun(out, AgeInt = AgeInt1)
	check       <- groupAges(newN, AgeS, AgeN = AgeN2)
	if (max(abs(check - Value2)) < tol){
		return(out)
	} else {
		rescaleAgeGroups(
				Value1 = out,  # only swap out Value1
				AgeInt1 = AgeInt1, 
				Value2 = Value2, 
				AgeInt2 = AgeInt2,
				splitfun = splitfun,
				tol = tol,
				recursive = recursive)
	}
}

#' force a (count) vector to abridged ages
#' @description This is a robustness utility, in place to avoid annoying hang-ups in \code{LTAbr()}. If data are given in non-standard ages, they are forced to standard abrdiged ages on the fly. Really this should happen prior to calling \code{LTAbr()}
#' @details This should be able to group up and group down as needed. \code{splitMono()} is used below the hood. \code{pclm()} or \code{splitUniform()} out to be flexible enough to do the same.
#' @inheritParams LTAbr
#' @seealso splitMono, LTAbr
#' @export
#' @examples
#' V1        <- pop1m_ind
#' Age       <- c(0,1,3,seq(5,100,5))
#' AgeInt    <- c(1,2,2,rep(5,19),1)
#' Value     <- tapply(V1,rep(Age,times=AgeInt), sum)
#' 
#' is_abridged(Age)
#' age_abridge_force(Value, AgeInt, Age)
age_abridge_force <- function(Value, AgeInt, Age){
  v1     <- splitMono(
    Value, 
    AgeInt,
    Age = Age)
  a1     <- min(Age):(length(v1)-1)
  AgeAbr <- calcAgeAbr(a1)
  vabr   <- tapply(v1, AgeAbr, sum)
  vabr
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




