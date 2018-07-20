
# Author: tim
###############################################################################

#' the basic Sprague age-splitting method
#' 
#' @description This method is based on the first stage of the Sprague R 
#' script prepared by Thomas Buettner and Patrick Gerland, itself based on the description
#' in Siegel and Swanson, 2004, p. 727.
#' 
#' @param popmat a numeric matrix of population counts in 5-year age groups, with integer-labeled 
#' margins (age in rows and year in columns).
#' @param Age integer vector. Lower age bound of age groups. Detected from row names of \code{popmat} if missing.
#' @param OAG logical (default \code{TRUE}. Is the final age group open?
#' @param closeout logical or character (default \code{"mono"}). 
#' @param pivotAge integer (default 90).
#' @param ... extra arguments passed to the closeout function.
#' @details Ages should refer to lower age bounds, ending in the open age group in the last row (not a closed terminal age). 
#' Dimension labelling is necessary. There must be at least six age groups (including the open group). One year of data will 
#' work as well, as long as it's given as a single-column matrix. There are different ways 
#' to specify closing out the graduation: To leave the Sprague results as-is specify \code{closeout = FALSE}. \code{TRUE} 
#' or \code{"mono"} will call the \code{monoCloseout()} function. The guarantees no negative values. Other closeout
#'  methods may be integrated in the future.
#' 
#' @return an age-period matrix od split population counts with the same number of 
#' columns as \code{popmat}, and single ages in rows.
#' 
#' @references 
#' \insertRef{sprague1880explanation}{DemoTools}
#' \insertRef{shryock1973methods}{DemoTools}
#' \insertRef{siegel2004methods}{DemoTools}
#' @export
#' 
#' @examples 
#' p5 <- structure(c(54170, 44775, 42142, 38464, 34406, 30386, 26933, 
#' 				23481, 20602, 16489, 14248, 9928, 8490, 4801, 3599, 2048, 941, 
#' 				326, 80, 17, 0, 57424, 44475, 41752, 39628, 34757, 30605, 27183, 
#' 				23792, 20724, 17056, 14059, 10585, 8103, 5306, 3367, 2040, 963, 
#' 				315, 80, 16, 1, 60272, 44780, 41804, 40229, 35155, 30978, 27456, 
#' 				24097, 20873, 17546, 13990, 11146, 7841, 5738, 3184, 2062, 961, 
#' 				311, 80, 15, 1, 62727, 45681, 42101, 40474, 35599, 31439, 27758, 
#' 				24396, 21055, 17958, 14046, 11589, 7731, 6060, 3086, 2083, 949, 
#' 				312, 79, 14, 1, 64816, 47137, 42508, 40532, 36083, 31940, 28092, 
#' 				24693, 21274, 18299, 14223, 11906, 7785, 6255, 3090, 2084, 938, 
#' 				316, 80, 14, 2), 
#' 		.Dim = c(21L, 5L), 
#' 		.Dimnames = list(seq(0,100,by=5), 1950:1954))
#' head(p5) # this is the entire matrix
#' # the last value is an open age group, preserve as such:
#' p1 <- sprague(p5, OAG = TRUE, closeout = FALSE)
#' head(p1); tail(p1)
#' colSums(p1) - colSums(p5) 
#' 
#' # another case, starting with single ages
#' # note sprague() does not group ages. You need to do it 
#' # first.
#' Value <- c(9544406,7471790,11590109,11881844,11872503,12968350,11993151,10033918,
#' 		14312222,8111523,15311047,6861510,13305117,7454575,9015381,10325432,
#' 		9055588,5519173,12546779,4784102,13365429,4630254,9595545,4727963,
#' 		5195032,15061479,5467392,4011539,8033850,1972327,17396266,1647397,
#' 		6539557,2233521,2101024,16768198,3211834,1923169,4472854,
#' 		1182245,15874081,1017752,3673865,1247304,1029243,12619050,1499847,
#' 		1250321,2862148,723195,12396632,733501,2186678,777379,810700,
#' 		7298270,1116032,650402,1465209,411834,9478824,429296,1190060,
#' 		446290,362767,4998209,388753,334629,593906,178133,
#' 		4560342,179460,481230,159087,155831,1606147,166763,93569,182238,
#' 		53567,1715697,127486,150782,52332,48664,456387,46978,34448,
#' 		44015,19172,329149,48004,28574,9200,7003,75195,13140,5889,
#' 		18915,21221,72373)
#' Age         <- 0:100
#' # group ages
#' Val5        <- groupAges(Value, Age)
#' # notice how this particular case produces a negative value in the last age
#' # before OAG:
#' pops <- sprague(popmat = Val5, OAG = TRUE, closeout = FALSE)
#' # this replaces ages 90+, guaranteed no negatives.
#' pops1 <- monoCloseout(popmat = Val5, pops = pops, OAG = TRUE)
#' # identical to:
#' pops2 <- sprague(popmat = Val5, OAG = TRUE, closeout = "mono")
#' stopifnot(all(pops1==pops2))

sprague <- function(popmat, Age, OAG = TRUE, closeout = "mono", pivotAge = 90){
	popmat            <- as.matrix(popmat)
	
	if (missing(Age)){
		Age               <- as.integer(rownames(popmat))
	}
	# this is innocuous if ages are already grouped
	pop5              <- apply(popmat, 2, groupAges, Age = Age, N = 5, shiftdown = 0)
	
	# generate coefficient matrix
	scm               <- spragueExpand(pop5, OAG = OAG)
	
	# redistribute
	pop1              <- scm %*% pop5
	
	# label and return
	zero              <- min(as.integer(rownames(popmat)))
	ages              <- zero:(nrow(scm)-1 + zero)
	dimnames(pop1)    <- list(ages, colnames(popmat))
	
	# default closeout with monoCloseout().
	# set to FALSE to turn off, write "mono"
	if (is.logical(closeout)){
		if (!closeout){
			return(pop1)
		}
		closeout <- "mono"
	}
	if (closeout == "mono"){
		pop1 <- monoCloseout(
				popmat = popmat, 
				pops = pop1, 
				OAG = OAG, 
				pivotAge = pivotAge)
	}
	pop1
}

#' create the Sprague coefficient matrix 
#' 
#' @description The resulting coefficient matrix is based on the number of rows in \code{popmat}
#' where we assume that each row of data is a 5-year age group. The final row may be an open 
#' or closed age group, as indicated by the \code{OAG} argument.
#' 
#' @param popmat numeric matrix of age-period population counts in 5-year age groups
#' @param OAG logical (default \code{TRUE}. Is the final age group open?
#' 
#' @details The \code{popmat} matrix is really just a placeholder in this case. This function is 
#' a utility called by the Sprague family of functions, where it is most convenient to just pass
#' in the same matrix being used in those calcs to determine the layout of the coefficient matrix.
#' 
#' @export
#' 
#' @references 
#' \insertRef{sprague1880explanation}{DemoTools}
#' \insertRef{shryock1973methods}{DemoTools}
#' \insertRef{siegel2004methods}{DemoTools}
#' @examples
#' p5 <- structure(c(54170, 44775, 42142, 38464, 34406, 30386, 26933, 
#' 				23481, 20602, 16489, 14248, 9928, 8490, 4801, 3599, 2048, 941, 
#' 				326, 80, 17, 0, 57424, 44475, 41752, 39628, 34757, 30605, 27183, 
#' 				23792, 20724, 17056, 14059, 10585, 8103, 5306, 3367, 2040, 963, 
#' 				315, 80, 16, 1, 60272, 44780, 41804, 40229, 35155, 30978, 27456, 
#' 				24097, 20873, 17546, 13990, 11146, 7841, 5738, 3184, 2062, 961, 
#' 				311, 80, 15, 1, 62727, 45681, 42101, 40474, 35599, 31439, 27758, 
#' 				24396, 21055, 17958, 14046, 11589, 7731, 6060, 3086, 2083, 949, 
#' 				312, 79, 14, 1, 64816, 47137, 42508, 40532, 36083, 31940, 28092, 
#' 				24693, 21274, 18299, 14223, 11906, 7785, 6255, 3090, 2084, 938, 
#' 				316, 80, 14, 2), 
#' 		.Dim = c(21L, 5L), 
#' 		.Dimnames = list(seq(0,100,by=5), 1950:1954))
#' coefsOA     <- spragueExpand(p5, TRUE)
#' coefsclosed <- spragueExpand(p5, FALSE)
#' dim(coefsOA)
#' dim(coefsclosed)

spragueExpand <- function(popmat, OAG = TRUE){
	popmat <- as.matrix(popmat)
	
	# figure out ages and years
	Age5   <- as.integer(rownames(popmat))
	Age1   <- min(Age5):max(Age5)
	yrs    <- as.integer(colnames(popmat))
	
	# nr 5-year age groups
	m      <- nrow(popmat)
	# nr rows in coef mat.
	n      <- m * 5 - ifelse(OAG, 4, 0)
	# number of middle blocks
	MP     <- m - ifelse(OAG, 5, 4) 
	
	# get the split coefficients
	# block for ages 0-9
	g1g2 <- matrix(c(
					0.3616, -0.2768,  0.1488, -0.0336,  0.0000, 
					0.2640, -0.0960,  0.0400, -0.0080,  0.0000, 
					0.1840,  0.0400, -0.0320,  0.0080,  0.0000, 
					0.1200,  0.1360, -0.0720,  0.0160,  0.0000, 
					0.0704,  0.1968, -0.0848,  0.0176,  0.0000, 
					0.0336,  0.2272, -0.0752,  0.0144,  0.0000, 
					0.0080,  0.2320, -0.0480,  0.0080,  0.0000, 
					-0.0080,  0.2160, -0.0080,  0.0000,  0.0000, 
					-0.0160,  0.1840,  0.0400, -0.0080,  0.0000, 
					-0.0176,  0.1408,  0.0912, -0.0144,  0.0000), 
			nrow = 10, ncol = 5, byrow = TRUE)
	# block for middle ages
	
	
	g3 <- matrix(c(-0.0128,   0.0848,  0.1504,   -0.0240,  0.0016, 
					-0.0016,   0.0144,  0.2224,   -0.0416,  0.0064, 
					0.0064,  -0.0336,  0.2544,   -0.0336,  0.0064, 
					0.0064,  -0.0416,  0.2224,    0.0144, -0.0016, 
					0.0016,  -0.0240,  0.1504,    0.0848, -0.0128),
			5, 5, byrow = TRUE) 
	
	# block prior to closeout
	g4g5 <- matrix(c(0.0000, -0.0144,  0.0912,  0.1408, -0.0176, 
					0.0000, -0.0080,  0.0400,  0.1840, -0.0160, 
					0.0000,  0.0000, -0.0080,  0.2160, -0.0080, 
					0.0000,  0.0080, -0.0480,  0.2320,  0.0080, 
					0.0000,  0.0144, -0.0752,  0.2272,  0.0336, 
					0.0000,  0.0176, -0.0848,  0.1968,  0.0704, 
					0.0000,  0.0160, -0.0720,  0.1360,  0.1200, 
					0.0000,  0.0080, -0.0320,  0.0400,  0.1840, 
					0.0000, -0.0080,  0.0400, -0.0960,  0.2640, 
					0.0000, -0.0336,  0.1488, -0.2768,  0.3616), 
			nrow = 10, ncol = 5, byrow = TRUE)
	
	
	
	## create a Sprague coefficient matrix for 5-year age groups
	bm               <- matrix(0, nrow = n, ncol =  m)
	## insert upper left block
	bm[1:10, 1:5]    <- g1g2
	
	# determine positions of middle blocks
	rowpos           <- matrix(11:((MP * 5) + 10), ncol = 5, byrow = TRUE)
	colpos           <- row(rowpos) + col(rowpos) - 1
	for (i in (1:MP)) {
		# calculate the slices and add middle panels accordingly
		bm[rowpos[i, ], colpos[i, ]] <- g3
	}
	
	## insert last two panels
	
	fr                <- nrow(bm) - ifelse(OAG,10,9)
	lr                <- fr + 9
	fc                <- ncol(bm) - ifelse(OAG, 5, 4)
	lc                <- fc + 4
	bm[fr:lr,fc:lc]   <- g4g5
	
	if (OAG){
		# preserve open ended age group
		bm[nrow(bm), ncol(bm)]    <- 1
	}
	
	bm
}


#' an oscillatory average of age splits
#' @description Single ages can be grouped into 5-year age groups in 5 ways by staggering terminal digits.
#' This method is a bit smoother than the standard Sprague or Beers methods, but not as smooth as \code{grabill()}.
#' 
#' @details This function works on a single vector of single-age counts, not on a matrix. Results are not
#' constrained to any particular age group, but are constrained to the total count. Negatives, \code{NA}, or \code{NaN} 
#' values are ignored in averaging. This can happen in older ages . It is recommended to run \code{monoCloseout()} or 
#' similar after the oscillatory split in such situations.
#' 
#' @param Value numeric vector of single age counts
#' @param Age integer vector of single ages (lower bound)
#' @param OAG logical (default \code{TRUE}). Is the last value the open age group?
#' @param splitfun function used to split at each digit grouping (default \code{sprague()}.
#' @param closeout logical or \code{"mono"} 
#' @param pivotAge Age to start blending in closeout values
#' @param ... optional arguments passed to \code{splitfun()}.
#' 
#' @return numeric vector of Sprague-smoothed counts
#' @references 
#' \insertRef{booth2015demographic}{DemoTools}
#' @export
#' @examples
#' # code currently breaking, needs to be revisited and updates completed, sorry
#' \dontrun{
#' Value <- structure(c(9544406, 7471790, 11590109, 11881844, 11872503, 12968350, 
#' 		 11993151, 10033918, 14312222, 8111523, 15311047, 6861510, 13305117, 
#' 		 7454575, 9015381, 10325432, 9055588, 5519173, 12546779, 4784102, 
#' 		 13365429, 4630254, 9595545, 4727963, 5195032, 15061479, 5467392, 
#' 		 4011539, 8033850, 1972327, 17396266, 1647397, 6539557, 2233521, 
#' 		 2101024, 16768198, 3211834, 1923169, 4472854, 1182245, 15874081, 
#' 		 1017752, 3673865, 1247304, 1029243, 12619050, 1499847, 1250321, 
#' 		 2862148, 723195, 12396632, 733501, 2186678, 777379, 810700, 7298270, 
#' 		 1116032, 650402, 1465209, 411834, 9478824, 429296, 1190060, 446290, 
#' 		 362767, 4998209, 388753, 334629, 593906, 178133, 4560342, 179460, 
#' 		 481230, 159087, 155831, 1606147, 166763, 93569, 182238, 53567, 
#' 		 1715697, 127486, 150782, 52332, 48664, 456387, 46978, 34448, 
#' 		 44015, 19172, 329149, 48004, 28574, 9200, 7003, 75195, 13140, 
#' 		 5889, 18915, 21221, 72373), .Names = 0:100)
#' #barplot(Value, main = "yup, these have heaping!")
#' # this is the basic case we compare with:
#' pop0    <- sprague(groupAges(Value),  OAG = TRUE)
#' # note: this function needs single ages to work because
#' # ages are grouped into 5-year age groups in 5 different ways.
#' # breaks
#' #pop1    <- splitOscillate(Value, OAG = TRUE, splitfun = sprague)
#' pop2    <- splitOscillate(Value, OAG = TRUE, splitfun = beers)
#' # what's smoother, splitOscillate() or grabill()?
#' # note, same closeout problem, can be handled by monoCloseout()
#' pop3    <- grabill(Value, OAG = TRUE)
#' # and technically you could give grabill as splitfun too
#' pop4   <- splitOscillate(Value, OAG = TRUE, splitfun = grabill)
#' 
#' Age <- 0:100
#' plot(Age, Value)
#' lines(Age, pop0, col = "blue")
#' # slightly smoother (also shifted though)
#' lines(Age, pop1)
#' # only different at very high ages, small nrs
#' lines(Age, pop2, col = "red", lty = 2, lwd = 2) 
#' lines(Age, pop3, col = "magenta")
#' lines(Age, pop4, col = "orange", lty = 2)
#' legend("topright", 
#' lty = c(1,1,2,1,2), 
#' lwd = c(1,1,2,1,1), 
#' col = c("blue","black","red","magenta","orange"),
#' 		legend = c("sprague()",
#'                 "splitOscillate(splitfun = sprague)", 
#' 				   "splitOscillate(splitfun = beers)",
#' 				   "grabill()",
#'                 "splitOscillate(splitfun = grabill)"))
#' 
#' # index of dissimilarity
#' ID(Value, pop0) # original vs sprague
#' ID(pop0,pop1) # sprague vs sprague osc
#' ID(pop1,pop2) # sprague osc vs beers osc
#' ID(pop2,pop3) # beers osc vs grabill
#' ID(pop3,pop4) # grabill vs grabill osc
#' # measre of smoothness:
#' mean(abs(diff(Value)))
#' mean(abs(diff(pop0)))
#' mean(abs(diff(pop1)))
#' mean(abs(diff(pop2)))
#' mean(abs(diff(pop3)))
#' mean(abs(diff(pop4)))
#' }
splitOscillate <- function(
		Value, 
		Age = 1:length(Value) - 1, 
		OAG = TRUE, 
		splitfun = sprague, 
		closeout = "mono", 
		pivotAge = 90, ...){
	
	N     <- length(Value)
	if (OAG){
		open   <- Value[N]
		OA     <- Age[N]
		Value  <- Value[-N]
		Age    <- Age[-N]
		N      <- N - 1
	} 
	TOT    <- sum(Value)
# select which ages to keep:
	p1x1   <- matrix(nrow = N, ncol = 5)
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
	
		# get first run estimate
		pop.est             <- splitfun(Val.i.5, OAG = FALSE, ...)
#        a                   <- rownames(pop.est)
#		if (closeout){
#			a.fake  <- (1:nrow(pop.est) - 1) * 5
#			pop.est <- monoCloseout(Val.i.5, Age = a.fake, pops = pop.est, OAG = FALSE)
#		}
		
		pop.est[pop.est < 0] <- 0
		p1x1[keep.i, i + 1]  <- pop.est
	}
	# take average per age
	p.out <- rowMeans(p1x1, na.rm = TRUE)
	# rescale to proper total
	p.out <- rescale.vector(p.out, TOT)
	# re-append the open age group if needed
	if (OAG){
		Age          <- c(Age, OA)
		p.out        <- c(p.out, open)
		names(p.out) <- Age
	}
	if (is.logical(closeout)){
		if (!closeout){
			return(p.out)
		}
		closeout <- "mono"
	}
	if (closeout == "mono"){
		p.out <- monoCloseout(popmat = Value, pops = p.out, OAG = OAG, pivotAge = 90)
	}
	
	p.out
}
