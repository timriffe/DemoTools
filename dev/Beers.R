

###############################################################################
# Based on BeersMSplit.R, whose header was:
########################################################################
## Interpolation of of event data (five year periods) 
## fifth/single years by Beers Modified six-term formula
## (Siegel and Swanson, 2004, p. 729)
## This formula applies some smoothing to the interpolant, which is
## recommended for time series of events with some dynamic
## R implementation by Thomas Buettner (21 Oct. 2015)
########################################################################
#ifn  <- "/home/tim/git/DemoTools/dev/Data/IND_Pop2015B1_1.7.csv"
#tp   <- t(read.csv(ifn, header = TRUE, row.names = 1, check.names = FALSE))
# popmat <- round(tp[,1:5])

#' create the Beers modified coefficient matrix 
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
#' @references 
#' \insertRef{beers1945modified}{DemoTools}
#' \insertRef{siegel2004methods}{DemoTools}
#' @export
#' popmat <- structure(c(54170, 44775, 42142, 38464, 34406, 30386, 26933, 
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
#' 
#' coefsOA     <- beersModExpand(popmat, OAG = TRUE)
#' coefsclosed <- beersModExpand(popmat, OAG = FALSE)
#' 
#' dim(beersModExpand(popmat, TRUE))
#' dim(beersModExpand(popmat, FALSE))
beersModExpand <- function(popmat, OAG = FALSE){
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
	
	## Beers Modified Split
	g1g2 <- matrix(c(
					 0.3332, -0.1938,  0.0702, -0.0118,  0.0022 ,
					 0.2569, -0.0753,  0.0205, -0.0027,  0.0006 ,
					 0.1903,  0.0216, -0.0146,  0.0032, -0.0005 ,
					 0.1334,  0.0969, -0.0351,  0.0059, -0.0011 ,
					 0.0862,  0.1506, -0.0410,  0.0054, -0.0012 ,
					
					 0.0486,  0.1831, -0.0329,  0.0021, -0.0009 ,
					 0.0203,  0.1955, -0.0123, -0.0031, -0.0004 ,
					 0.0008,  0.1893,  0.0193, -0.0097,  0.0003 ,
					-0.0108,  0.1677,  0.0577, -0.0153,  0.0007 ,
					-0.0159,  0.1354,  0.0972, -0.0170,  0.0003
			), nrow = 10, ncol = 5, byrow = TRUE)
	
	g3 <- matrix(c(
					-0.0160,  0.0973,  0.1321, -0.0121, -0.0013 ,
					-0.0129,  0.0590,  0.1564,  0.0018, -0.0043 ,
					-0.0085,  0.0260,  0.1650,  0.0260, -0.0085 ,
					-0.0043,  0.0018,  0.1564,  0.0590, -0.0129 ,
					-0.0013, -0.0121,  0.1321,  0.0973, -0.0160
			), 5, 5, byrow = TRUE) 
	
	g4g5 <- matrix(c(0.0003, -0.0170,  0.0972,  0.1354, -0.0159 ,
					 0.0007, -0.0153,  0.0577,  0.1677, -0.0108 ,
					 0.0003, -0.0097,  0.0193,  0.1893,  0.0008 ,
					-0.0004, -0.0031, -0.0123,  0.1955,  0.0203 ,
					-0.0009,  0.0021, -0.0329,  0.1831,  0.0486 ,
					
					-0.0012,  0.0054, -0.0410,  0.1506,  0.0862 ,
					-0.0011,  0.0059, -0.0351,  0.0969,  0.1334 ,
					-0.0005,  0.0032, -0.0146,  0.0216,  0.1903 ,
					 0.0006, -0.0027,  0.0205, -0.0753,  0.2569 ,
					 0.0022, -0.0118,  0.0702, -0.1938,  0.3332 
			), nrow = 10, ncol = 5, byrow = TRUE)
	
	
	## create a Beers coefficient matrix for 5-year age groups
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
	## standard format for Beers coefficients
	
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