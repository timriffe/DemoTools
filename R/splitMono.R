
#' split age groups using a monotonic spline
#' @description Take the cumulative sum of \code{Value} and then run a monotonic spline through it. The first 
#' differences split back single-age estimates of \code{Value}. Optionally keep the open age group untouched. 
#' 
#' @details We use the \code{"monoH.FC"} method of \code{stats::splinefun()} to fit the spline because 1)
#' it passes exactly through the points, 2) it is monotonic and therefore guarantees positive counts, and 3) 
#' it seems to be a bit less wiggly (lower average first differences of split counts) than a pchip tends to do, 
#' at least in the tested data. Data may be given in single or grouped ages.
#' 
#' @param Value numeric vector of counts in age groups
#' @param Age integer vector of lower bound of age groups
#' @param OAG logical (default \code{FALSE}). Would we like to re-impute the last 
#' element of \code{Value} as the open age group?
#' @return numeric vector of single age counts 
#' @importFrom stats splinefun
#' @references 
#' \insertRef{fritsch1980monotone}{DemoTools}
#' @export
#' @examples
#' Value <- structure(c(88623, 90842, 93439, 96325, 99281, 102051, 104351, 
#'				 106555, 109170, 112188, 113582, 112614, 108904, 102622, 95867, 
#'				 80874, 60196, 37523, 17927, 5642, 1110), .Names = c("0", "5", 
#'				 "10", "15", "20", "25", "30", "35", "40", "45", "50", "55", "60", 
#'				 "65", "70", "75", "80", "85", "90", "95", "100"))
#'
#' # overwrite open age group with a single age estimate for that age
#' # (doesn't extrapolate)
#' splitMono(Value)
#' # or respect open age group
#' splitMono(Value, OAG = TRUE)
#' 
#' # Also accepts single ages:
#' Value <- structure(c(9544406, 7471790, 11590109, 11881844, 11872503, 12968350, 
#' 		11993151, 10033918, 14312222, 8111523, 15311047, 6861510, 13305117, 
#' 		7454575, 9015381, 10325432, 9055588, 5519173, 12546779, 4784102, 
#' 		13365429, 4630254, 9595545, 4727963, 5195032, 15061479, 5467392, 
#' 		4011539, 8033850, 1972327, 17396266, 1647397, 6539557, 2233521, 
#' 		2101024, 16768198, 3211834, 1923169, 4472854, 1182245, 15874081, 
#' 		1017752, 3673865, 1247304, 1029243, 12619050, 1499847, 1250321, 
#' 		2862148, 723195, 12396632, 733501, 2186678, 777379, 810700, 7298270, 
#' 		1116032, 650402, 1465209, 411834, 9478824, 429296, 1190060, 446290, 
#' 		362767, 4998209, 388753, 334629, 593906, 178133, 4560342, 179460, 
#' 		481230, 159087, 155831, 1606147, 166763, 93569, 182238, 53567, 
#' 		1715697, 127486, 150782, 52332, 48664, 456387, 46978, 34448, 
#' 		44015, 19172, 329149, 48004, 28574, 9200, 7003, 75195, 13140, 
#' 		5889, 18915, 21221, 72373), .Names = 0:100)
#'  
#'  from_single  <- splitMono(Value, OAG = TRUE) 
#'  # compare, had we pre-grouped ages:
#'  Val5         <- groupAges(Value, N=5, shiftdown = 0)
#'  from_grouped <- splitMono(Val5, OAG = TRUE)
#'  # de facto unit test:
#'  stopifnot(all(from_single == from_grouped))

splitMono <- function(Value, Age, OAG = FALSE){
	if (missing(Age)){
		Age               <- as.integer(names(Value))
	}
	
	# this is innocuous if ages are already grouped
	Age5       <- calcAgeN(Age, N = 5)
	pop5       <- groupAges(Value, AgeN = Age5, shiftdown = 0)
	
	AgePred    <- min(Age):(max(Age) + 1)
	y          <- c(0, cumsum(pop5))
	x          <- c(0, sort(unique(Age5)) + 5)
	y1         <- splinefun(y ~ x, method = "monoH.FC")(AgePred)
	single.out <- diff(y1)
	if (OAG){
		single.out[length(single.out)] <- pop5[length(pop5)]
	}
	single.out
}

#' blend the Sprague upper boundary age estimates into monotonic spline estimates
#' 
#' @description A simple monotonic spline on the cumulative sum of population counts
#' may return more convincing single age count estimates than the Sprague or other splitting methods.
#' This function blends the given single age population estimates starting at \code{pivotAge}.
#' 
#' @param popmat a numeric matrix of population counts in 5-year age groups, with integer-labeled 
#' margins (age in rows and year in columns).
#' @param pops optional numeric matrix of single age population counts derived from \code{popmat}.
#' @param pivotAge integer (default 90). Age at which to switch to spline-based estimates.
#' @param splitfun optional. The function used to create pops. Default \code{sprague}. 
#' Could also be \code{grabill}, \code{beersModSimple}, or any other function that similarly transforms.
#' @param OAG logical (default \code{FALSE}). Would we like to re-impute the last 
#' element of \code{Value} as the open age group?
#' @param ... arguments to be optionally passed to \code{splitfun()}.
#' @return numeric matrix of age by year estimates of single-age counts.
#' 
#' @details The \code{pivotAge} must be at least 10 years below the maximum age detected from 
#' \code{rownames(popmat)}, but not lower than 80. In the exact \code{pivotAge}, we may either take
#' the Sprague estimates or the spline estimates, depending on which is larger, then the single-age estimates
#' for this 5-year age group are rescaled to sum to the original total in \code{popmat}. Higher ages are taken from
#' the spline-based age splits. The spline results are derive from the \code{"monoH.FC"} method of \code{splinefun()} 
#' on the cumulative sum of the original age grouped data. One could use this function to perform the same
#' closeout to Grabill estimates, if these are given via the \code{pops} argument. See examples. Note
#' that the Grabill split method mixed with this closeout will not necessarily preserve the annual totals,
#' and this function performs to rescaling. The open age group is preserved (and must be included in \code{popmat}).
#' 
#' @export 
#' 
#' @examples
#'  popmat <- structure(c(54170, 44775, 42142, 38464, 34406, 30386, 26933, 
#' 23481, 20602, 16489, 14248, 9928, 8490, 4801, 3599, 2048, 941, 
#' 326, 80, 17, 0, 57424, 44475, 41752, 39628, 34757, 30605, 27183, 
#' 23792, 20724, 17056, 14059, 10585, 8103, 5306, 3367, 2040, 963, 
#' 315, 80, 16, 1, 60272, 44780, 41804, 40229, 35155, 30978, 27456, 
#' 24097, 20873, 17546, 13990, 11146, 7841, 5738, 3184, 2062, 961, 
#' 311, 80, 15, 1, 62727, 45681, 42101, 40474, 35599, 31439, 27758, 
#' 24396, 21055, 17958, 14046, 11589, 7731, 6060, 3086, 2083, 949, 
#' 312, 79, 14, 1, 64816, 47137, 42508, 40532, 36083, 31940, 28092, 
#' 24693, 21274, 18299, 14223, 11906, 7785, 6255, 3090, 2084, 938, 
#' 316, 80, 14, 2), .Dim = c(21L, 5L), .Dimnames = list(c("0", "5", 
#' "10", "15", "20", "25", "30", "35", "40", "45", "50", "55", "60", 
#' "65", "70", "75", "80", "85", "90", "95", "100"), c("1950", "1951", 
#' "1952", "1953", "1954")))
#' 
#' #closed.out <- monoCloseout(popmat)
#' #colSums(closed.out) - colSums(popmat)
#' #monoCloseout(popmat, pivotAge = 85)
#' # giving a different single-age split to close out this way:
#' popg <- grabill(popmat)
#' grabill.closed.out <- monoCloseout(popmat, popg)
#' # totals not necessarily preserved if mixed w Grabill
#' # I wouldn't recommend a rescale of the total, since the 
#' # only part we mess with here is the old age section. Ergo,
#' # one may wish to instead rescale results colSums() of 
#' # popg at age pivotAge and higher.
#' colSums(grabill.closed.out) - colSums(popmat)
#' # also works on an age-labelled vector of data
#' popvec <- popmat[,1]
#' closed.vec <- monoCloseout(popvec)
#' # let's compare this one with sprague()
#' simple.vec <- sprague(popvec)
#' # and with a simple monotonic spline
#' mono.vec <- splitMono(popvec)
#' \dontrun{
#' plot(85:100,simple.vec[86:101], type = 'l', main = "In this case sprague() is the smoothest")
#' lines(85:100,closed.vec[86:101], col = "red", lwd = 2)
#' lines(85:100,mono.vec[86:101], col = "blue", lty = 2)
#' legend("topright",lty=c(1,2,2), col = c("black","red","blue"),lwd = c(1,2,1),
#' 		legend = c("sprague()","monoCloseout()", "splitMono()"))
#' }
# TODO: make this deal w OAG in a more consistent way
monoCloseout <- function(popmat, pops, pivotAge = 90, splitfun = sprague, OAG = FALSE, ...){
	popmat <- as.matrix(popmat)
	if (missing(pops)){
		pops    <- splitfun(popmat, ...)
	}
	# get the spline population split
	popmono <- apply(popmat, 2, splitMono, OAG = OAG)
	
	# some age pars
	Age5    <- as.integer(rownames(popmat))
	Age1    <- min(Age5):max(Age5)
	
	# some checks on pivotAge...
	if (!(max(Age1) - 10) >= pivotAge){
		pivotAge <- max(Age1) - 10
		if (pivotAge < 80){
			warning("pivotAge wasn't in rownames(popmat), moved it to 3rd 
							from bottom row of popmat, but appears to be < 80
							so returning sprague() output as-is, no extra closeout performed.")
			return(pops)
		}
		warning("pivotAge moved to ", pivotAge, ", continued.")
	}
	# -----------------------------
	# now begin the closeout blend.
	p.i              <- which(Age1 == pivotAge)
	## substitute Sprague interpolation if > pchip for better belnding of the two series
	pop.c            <- popmono[p.i:(p.i + 4), , drop = FALSE]
	ind              <- pops[p.i, ] > pop.c[1, ]
	pop.c[1, ind]    <- pops[p.i, ind] 
	
	## adjust back on initial pop 5x5 for age 90-94
	## proportional distribution
	pop.c[is.na(pop.c )] <- 0
	prop             <- prop.table(pop.c, margin = 2)
	pivot5           <- popmat[as.character(pivotAge), ]
	pop.c            <- t(t(prop) * pivot5)
	## append the remaining of the age groups (except last open age)
	## 95-99 onward
	m                <- nrow(pops)
	pop.c            <- rbind(pop.c, popmono[(p.i + 5):m, , drop = FALSE])
	## append Sprague interpolation before age 90
	pop.c            <- rbind(pops[1:(p.i - 1), , drop = FALSE], pop.c)
	
	## deal with negative values if applicable (but in principle should not be happening)
	pop.c[pop.c < 0] <- 0
	
	# label and return
	#dimnames(pop.c) <- list(Age1, colnames(popmat))
	rownames(pop.c) <- Age1
	colnames(pop.c) <- colnames(popmat)
	pop.c
}




