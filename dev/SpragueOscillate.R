
# Author: tim
###############################################################################

#' an oscillatory average of Sprague age splits
#' @description Single ages can be grouped into 5-year age groups in 5 ways by staggering terminal digits.
#' This method is a bit smoother than the standard \code{spragueSimple()} method, but not as smooth as \cide{grabill()}.
#' 
#' @details This function works on a single vector of single-age counts, not on a matrix. Results are not
#' constrained to any particular age group, but are constrained to the total count.
#' The option to closeout using \code{spragueCloseout()} is recommended because it usually gives 
#' more plausible results and it avoids negative values. This is run separately on each Sprague split,
#' rather than on the aggregate results. 
#' 
#' @param Value numeric vector of single age counts
#' @param Age integer vector of single ages (lower bound)
#' @param OAG logical (default \code{TRUE}). Is the last value the open age group?
#' @param closeout logical (default \code{TRUE}). Shall we close out each sprague split with a monotonic spline fit?
#' 
#' @return numeric vector of Sprague-smoothed counts
#' @export
#' @examples
#' Value <- c(9544406,7471790,11590109,11881844,11872503,12968350,11993151,10033918,
#' 14312222,8111523,15311047,6861510,13305117,7454575,9015381,10325432,
#' 9055588,5519173,12546779,4784102,13365429,4630254,9595545,4727963,
#' 5195032,15061479,5467392,4011539,8033850,1972327,17396266,1647397,
#' 6539557,2233521,2101024,16768198,3211834,1923169,4472854,
#' 1182245,15874081,1017752,3673865,1247304,1029243,12619050,1499847,
#' 1250321,2862148,723195,12396632,733501,2186678,777379,810700,
#' 7298270,1116032,650402,1465209,411834,9478824,429296,1190060,
#' 446290,362767,4998209,388753,334629,593906,178133,
#' 4560342,179460,481230,159087,155831,1606147,166763,93569,182238,
#' 53567,1715697,127486,150782,52332,48664,456387,46978,34448,
#' 44015,19172,329149,48004,28574,9200,7003,75195,13140,5889,
#' 18915,21221,72373)
#' Age <- 0:100
#' names(Value) <- Age
#' #barplot(Value, main = "yup, these have heaping!")
#' # this is the basic case we compare with:
#' pop0    <- spragueSimple(groupAges(Value,Age))
#' # note: this function needs single ages to work because
#' # ages are grouped into 5-year age groups in 5 different ways.
#' (pop1   <- spragueOscillate(Value, Age, closeout = FALSE))
#' # see the NaN value? That because there were some negatives produced by 
#' # spragueSimple(). We can call spragueCloseout() inside spragueOscillate()
#' # to handle such cases:
#' (pop2   <- spragueOscillate(Value, Age, closeout = TRUE))
#' # what's smoother, spragueOscillate() or grabill()?
#' # note, same closeout problem, can be handled by spragueCloseout()
#' (pop3   <- grabill(groupAges(Value, Age)))
#' #pop4   <- spragueCloseout(groupAges(Value, Age), pops = pop3)
#' \dontrun{
#' plot(Age, Value)
#' lines(Age, pop0, col = "blue")
#' # slightly smoother (also shifted though)
#' lines(Age, pop1)
#' # only different at very high ages, small nrs
#' lines(Age, pop2, col = "red", lty = 2, lwd = 2) 
#' lines(Age, pop3, col = "magenta")
#' legend("topright", lty = c(1,1,2,1), lwd = c(1,1,2,1), col = c("blue","black","red","magenta"),
#' 		legend = c("spragueSimple()",
#'                 "spragueOscillate(closeout = FALSE)", 
#' 				   "spragueOscillate(closeout = TRUE)",
#' 				   "grabill()"))
#' }

spragueOscillate <- function(Value, Age, OAG = TRUE, closeout = TRUE){
	
	N     <- length(Value)
	if (OAG){
		open   <- Value[N]
		OA     <- Age[N]
		Value  <- Value[-N]
		Age    <- Age[-N]
		N      <- N - 1
	} 
	TOT <- sum(Value)
# select which ages to keep:
	p1x1   <- matrix(nrow = length(Value), ncol = 5)
	rownames(p1x1) <- Age
	for (i in 0:4){
		# regroup ages
		Age.i.5             <- calcAgeN(Age, shiftdown = i)
		# only use age groups w 5 single ages represented
		keep.i              <- rep(rle(Age.i.5)$leng, rle(Age.i.5)$leng) == 5
		# cut vector down to those cases
		Age.i.5             <- Age.i.5[keep.i]
		# cut counts down to those cases
		Val.i               <- Value[keep.i]
		# group ages into said 5-year age groups
		Val.i.5             <- groupAges(Val.i, AgeN = Age.i.5)
		# make fake open age
		Val.i.5             <- c(Val.i.5, pi)
		names(Val.i.5)      <- c(unique(Age.i.5), max(Age.i.5) + 5)
		# get first run estimate
		pop.est             <- spragueSimple(Val.i.5)
		if (closeout){
			pop.est <- spragueCloseout(Val.i.5, pop.est)
		}
		
		pop.est[pop.est < 0] <- NA
		pop.est             <- pop.est[-length(pop.est)]
		p1x1[keep.i, i + 1] <- pop.est
	}
	# take average per age
	p.out <- rowMeans(p1x1, na.rm = TRUE)
	# rescale to proper total
	p.out <- p.out * TOT / sum(p.out, na.rm = TRUE)
	# re-append the open age group if needed
	if (OAG){
		Age   <- c(Age, OA)
		p.out <- c(p.out, open)
		names(p.out) <- Age
	}
	
	p.out
}


