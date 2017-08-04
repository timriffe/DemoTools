

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
# Author: tim
tp <- matrix(0, nrow = 15, ncol = 5)
## dimensions
NCOL <- dim(tp)[2]   ## countries or trajectories
NAG5 <- dim(tp)[1]   ## age or time

# number of single year or single age groups
NAG1 <- NAG5 * 5

beersModExpand <- function(popmat, OAG = FALSE){
	
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
	
	bm               <- matrix(0, nrow = m1 + 1, ncol =  m)
	
	## insert upper left block
	bm[1:10, 1:5]    <- g1g2
	
	# determine positions of middle blocks
	rowpos           <- matrix(11:((MP*5) + 10), ncol = 5, byrow = TRUE)
	colpos           <- row(rowpos) + col(rowpos) - 1
	for (i in (1:MP)) {
		# calculate the slices and add middle panels accordingly
	    bm[rowpos[i,], colpos[i, ]] <- g3
	}
	## standard format for Beers coefficients
	
	bm[,]
}