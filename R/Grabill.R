#' create the Grabill coefficient matrix 
#' 
#' @description The resulting coefficient matrix is based on the number of rows in \code{popmat}
#' where we assume that each row of data is a 5-year age group and the final row is an open age group
#' to be preserved as such.
#' 
#' @param popmat numeric matrix of age-period population counts in 5-year age groups
#' @param OAG logical (default \code{TRUE}. Is the final age group open?
#' 
#' @details The \code{popmat} matrix is really just a placeholder in this case. This function is 
#' a utility called by the Grabill family of functions, where it is most convenient to just pass
#' in the same matrix being used in those calcs to determine the layout of the coefficient matrix.
#' Note that these coefficients do not constrain population counts to their year totals. This function 
#' is called by \code{grabill()}, which ensures matching marginals by 1) blending boundary ages 
#' into the Sprague estimated population, and 2) a second constraint on the middle age groups to enforce
#' matching sums.
#' 
#' @references
#' \insertRef{shryock1973methods}{DemoTools}
#' 
#' @export
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
#' grabillExpand(p5, OAG = TRUE)
#' grabillExpand(p5, OAG = FALSE)
grabillExpand <- function(popmat, OAG = TRUE){
	popmat            <- as.matrix(popmat)
	
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
	
	# primary grabill coef block
	g3g              <- matrix(
			c(
					0.0111,	0.0816,	 0.0826,	0.0256,	-0.0009,
					0.0049,	0.0673,	 0.0903,	0.0377,	-0.0002,
					0.0015,	0.0519,	 0.0932,	0.0519,	 0.0015,
					-0.0002,	0.0377,	 0.0903,	0.0673,	 0.0049,
					-0.0009,	0.0256,	 0.0826,	0.0816,	 0.0111),
			5, 5, byrow = TRUE)
	## create a Grabill coefficient matrix for 5-year age groups
	gm               <- matrix(0, nrow = n, ncol =  m)
	
	
	fr                <- nrow(gm) - ifelse(OAG,10,9)
	lr                <- fr + 9
	fc                <- ncol(gm) - ifelse(OAG, 5, 4)
	lc                <- fc + 4
	
	# ----------------------------------------------------------
	# Note: for the boundary ages we keep shuffling in g3g, the same grabill
	# coefs. The columns on the boundaries will NOT sum to 1. These coefs are
	# used just for the firs pass, then results blend into the Sprague boundary
	# estimates.
	# ----------------------------------------------------------
	# the young age coefficients
	g1g2g              <- matrix(0,nrow=10,ncol=5)
	g1g2g[1:5, 1:3]    <- g3g[,3:5]
	g1g2g[6:10, 1:4]   <- g3g[,2:5]
	# the old age coefficients
	g4g5g              <- matrix(0,nrow=10,ncol=5)
	g4g5g[1:5, 2:5]    <- g3g[,1:4]
	g4g5g[6:10, 3:5]   <- g3g[,1:3]
	
	gm[1:10, 1:5]    <- g1g2g
	gm[fr:lr,fc:lc]  <- g4g5g
	
	
	# determine positions of middle blocks
	rowpos             <- matrix(11:((MP*5) + 10), ncol = 5, byrow = TRUE)
	colpos             <- row(rowpos) + col(rowpos) - 1
	for (i in (1:MP)) {
		# calculate the slices and add middle panels accordingly
		gm[rowpos[i,], colpos[i, ]] <- g3g
	}
	
	if (OAG){
		# preserve open ended age group
		gm[nrow(gm), ncol(gm)]    <- 1
	}
	
	# return coefficient matrix
	gm
}

#' the basic Grabill age-splitting method
#' 
#' @description This method uses Grabill's aggressive redistribution of middle ages and blends into
#' Sprague estimated single-age population counts for the first and final ten ages. Open age groups
#' are preserved, as are annual totals.
#' 
#' @param popmat a numeric matrix of population counts in 5-year age groups, with integer-labeled 
#' margins (age in rows and year in columns).
#' @param Age integer vector lower age bound of age groups
#' @param OAG logical (default \code{TRUE}. Is the final age group open?
#' @details Ages should refer to lower age bounds, ending in the open age group in the last row (not a closed terminal age). 
#' Dimension labelling is necessary. There must be at least six age groups (including the open group). One year of data will 
#' work as well, as long as it's given as a single-column matrix. Data may be given in either single or grouped ages.
#' 
#' @return an age-period matrix od split population counts with the same number of 
#' columns as \code{popmat}, and single ages in rows.
#' 
#' @references 
#' \insertRef{shryock1973methods}{DemoTools}
#' 
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
#' p1g <- grabill(p5)
#' head(p1g); tail(p1g)
#' colSums(p1g) - colSums(p5) 
#' p1s <- spragueSimple(p5)
#' \dontrun{
#' plot(seq(0,100,by=5),p5[,1]/5,type = "s", col = "gray", xlab = "Age", ylab = "Count")
#' lines(0:100, p1g[,1], col = "red", lwd = 2)
#' lines(0:100, p1s[,1], col = "blue", lty = 2)
#' legend("topright", 
#'		lty = c(1,1,2), 
#'		col = c("gray","red","blue"), 
#'		lwd = c(1,2,1), 
#'		legend = c("grouped","Grabill", "Sprague"))
#' }
#' 
#' # also works for single ages:
#'  popmat <- structure(
#'     c(9544406, 7471790, 11590109, 11881844, 11872503, 12968350, 
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
#' 		5889, 18915, 21221, 72373), .Dim = c(101L, 1L), .Dimnames = list(
#' 		0:100, NULL))
#' grab1 <- grabill(popmat)
#' \dontrun{
#' plot(0:100, c(popmat))
#' lines(0:100, c(grab1))
#' }
grabill <- function(popmat, Age, OAG = TRUE){
	popmat            <- as.matrix(popmat)
	
	if (missing(Age)){
		Age               <- as.integer(rownames(popmat))
	}
	# this is innocuous if ages are already grouped
	pop5              <- apply(popmat, 2, groupAges, Age = Age, N = 5, shiftdown = 0)
	
	
	# get coefficient matrices for Sprague and Grabill
	scmg              <- grabillExpand(pop5, OAG = OAG)
	scm               <- spragueExpand(pop5, OAG = OAG)
	
	# split pop counts
	pops              <- scm %*% pop5
	popg              <- scmg %*% pop5
	
	# ---------------------------------------------
	# now we graft the two estimates in together,
	# preserving the middle part for grabill, and blending
	# aggressively into the young and closeout parts of Sprague
	# weights for grafting in grabill
	m                 <- nrow(pops)
	lr                <- m - 1
	fr                <- lr - 9
	
	# these weights do much better than linear weights.
	w10               <- exp(row(pops[1:10, , drop = FALSE]) ) / exp(10.1)
	
	# blend together young ages
	popg[1:10, ]      <- w10 * popg[1:10, ] + (1 - w10) * pops[1:10, ]
	
	# blend together old ages
	popg[fr:lr, ]     <- w10[10:1, ] * popg[fr:lr, ] + (1 - w10[10:1, ]) * pops[fr:lr, ]
	
	# ---------------------------------------------
	# now we take care of the marginal constraint problem
	# make weighting matrix 
	wr                <- pops * 0 + 1
	wr[1:10, ]        <- w10
	wr[fr:lr, ]       <- w10[10:1, ]
	wr[nrow(wr), ]    <- 0
	
	# weighted marginal sums. The difference we need to redistribute
	redist            <- colSums(pops) - colSums(popg)
	
	middle.part       <- popg * wr
	
	# the difference to redistribute
	add.in            <- t(t(prop.table(middle.part,2)) * redist)
	popg              <- popg + add.in
	# ---------------------------------------------
	# label dims and return
	rg                <- range(as.integer(rownames(popmat)))
	dimnames(popg)    <- list(rg[1]:(rg[2]+ifelse(OAG, 0, 4)), colnames(popmat))
	
	popg
}



