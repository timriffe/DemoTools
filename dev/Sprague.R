
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
#' @details Ages should refer to lower age bounds, ending in the open age group in the last row (not a closed terminal age). 
#' Dimension labelling is necessary. There must be at least six age groups (including the open group). One year of data will 
#' work as well, as long as it's given as a single-column matrix.
#' 
#' @return an age-period matrix od split population counts with the same number of 
#' columns as \code{popmat}, and single ages in rows.
#' 
#' @references 
#' Shryock, H. S., Siegel, J. S., & Larmon, E. A. (1973). 
#' The methods and materials of demography. US Bureau of the Census.
#' 
#' Seigel, J. S., & Swanson, D. A. (2004). T
#' he methods and materials of demography. Elsevier Academic Press, London.
#' @export
#' 
#' @examples 
#' p5 <- structure(c(54170.08, 44774.6, 42141.587, 38463.515, 34405.607, 
#' 162369.816, 57424.3568738, 44475.4981681, 41751.7574114, 39628.4338929, 
#' 34756.9473002, 164194.0485702, 60272.2061248, 44780.1982856, 
#' 41803.6541424, 40229.0292664, 35154.7682192, 166275.9022992, 
#' 62726.896388, 45681.1355532, 42100.72506, 40473.8600572, 35598.545404, 
#' 168556.5331816, 64815.5458002, 47136.5341033, 42508.3026466, 
#' 40532.3096745, 36082.7490698, 170990.473735, 66579.122, 49070.407, 
#' 42953.604, 40534.586, 36596.844, 173545.633), .Dim = c(6L, 6L
#' ), .Dimnames = list(seq(0,25,5), 1950:1955))
#' head(p5) # this is the entire matrix
#' p1 <- spragueSimple(p5)
#' head(p1); tail(p1)
#' colSums(p1) - colSums(p5) 

spragueSimple <- function(popmat){
	popmat            <- as.matrix(popmat)
	scm               <- spragueExpand(popmat)
	
	pop1              <- scm %*% popmat
	
	rg <- range(as.integer(rownames(popmat)))
	dimnames(pop1)    <- list(rg[1]:rg[2], colnames(popmat))
	pop1
}

#' create the Sprague coefficient matrix 
#' 
#' @description The resulting coefficient matrix is based on the number of rows in \code{popmat}
#' where we assume that each row of data is a 5-year age group and the final row is an open age group
#' to be preserved as such.
#' 
#' @param popmat numeric matrix of age-period population counts in 5-year age groups
#' 
#' @details The \code{popmat} matrix is really just a placeholder in this case. This function is 
#' a utility called by the Sprague family of functions, where it is most convenient to just pass
#' in the same matrix being used in those calcs to determine the layout of the coefficient matrix.
#' 
#' @export

spragueExpand <- function(popmat){
	popmat <- as.matrix(popmat)
	
	# figure out ages and years
	Age5   <- as.integer(rownames(popmat))
	Age1   <- min(Age5):max(Age5)
	yrs    <- as.integer(colnames(popmat))
	
	# nr 5-year age groups
	m      <- nrow(popmat)
	# nr closed single ages
	m1     <- m * 5 - 5 
	# number of middle blocks
	MP     <- m - 5 
	
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
	scm               <- matrix(0, nrow = m1 + 1, ncol =  m)
	
	## insert upper left block
	scm[1:10, 1:5]    <- g1g2
	
	# determine positions of middle blocks
	rowpos           <- matrix(11:((MP*5) + 10), ncol = 5, byrow = TRUE)
	colpos           <- row(rowpos) + col(rowpos) - 1
	for (i in (1:MP)) {
		# calculate the slices and add middle panels accordingly
		scm[rowpos[i,], colpos[i, ]] <- g3
	}
	
	## insert last two panels
	fr                <- (m - 3) * 5 + 1
	lr                <- (m - 1) * 5
	fc                <- MP 
	lc                <- MP + 4 
	scm[fr:lr,fc:lc]  <- g4g5
	
	# last open ended age group
	scm[m1 + 1, m]    <- 1
	
	scm
}


#' create the Grabill coefficient matrix 
#' 
#' @description The resulting coefficient matrix is based on the number of rows in \code{popmat}
#' where we assume that each row of data is a 5-year age group and the final row is an open age group
#' to be preserved as such.
#' 
#' @param popmat numeric matrix of age-period population counts in 5-year age groups
#' 
#' @details The \code{popmat} matrix is really just a placeholder in this case. This function is 
#' a utility called by the Grabill family of functions, where it is most convenient to just pass
#' in the same matrix being used in those calcs to determine the layout of the coefficient matrix.
#' Note that these coefficients do not constrain population counts to their year totals. This function 
#' is called by \code{grabill()}, which ensures matching marginals by 1) blending boundary ages 
#' into the Sprague estimated population, and 2) a second constraint on the middle age groups to enforce
#' matching sums.
#' 
#' @export

grabillExpand <- function(popmat){
	popmat            <- as.matrix(popmat)

	# nr 5-year age groups
	m                 <- nrow(popmat)
	# nr closed single ages
	m1                <- m * 5 - 5 
	# number of middle blocks
	MP                <- m - 5 
	
	## create a Grabill coefficient matrix for 5-year age groups
	scmg              <- matrix(0, nrow = m1 + 1, ncol =  m)
	
	## insert last two panels
	fr                <- (m - 3) * 5 + 1
	lr                <- (m - 1) * 5
	fc                <- MP 
	lc                <- MP + 4 
	
	# preserve open age group
	scmg[m1 + 1, m]   <- 1
	
	# primary grabill coef block
	g3g              <- matrix(
			             c(
					       0.0111,	0.0816,	 0.0826,	0.0256,	-0.0009,
					       0.0049,	0.0673,	 0.0903,	0.0377,	-0.0002,
					       0.0015,	0.0519,	 0.0932,	0.0519,	 0.0015,
					      -0.0002,	0.0377,	 0.0903,	0.0673,	 0.0049,
					      -0.0009,	0.0256,	 0.0826,	0.0816,	 0.0111),
			             5, 5, byrow = TRUE)
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
	
	scmg[1:10, 1:5]    <- g1g2g
	scmg[fr:lr,fc:lc]  <- g4g5g
	
	# determine positions of middle blocks
	rowpos             <- matrix(11:((MP*5) + 10), ncol = 5, byrow = TRUE)
	colpos             <- row(rowpos) + col(rowpos) - 1
	for (i in (1:MP)) {
		# calculate the slices and add middle panels accordingly
		scmg[rowpos[i,], colpos[i, ]] <- g3g
	}
	# return coefficient matrix
	scmg
}

#' the basic Grabill age-splitting method
#' 
#' @description This method uses Grabill's aggressive redistribution of middle ages and blends into
#' Sprague estimated single-age population counts for the first and final ten ages. Open age groups
#' are preserved, as are annual totals.
#' 
#' @param popmat a numeric matrix of population counts in 5-year age groups, with integer-labeled 
#' margins (age in rows and year in columns).
#' @details Ages should refer to lower age bounds, ending in the open age group in the last row (not a closed terminal age). 
#' Dimension labelling is necessary. There must be at least six age groups (including the open group). One year of data will 
#' work as well, as long as it's given as a single-column matrix.
#' 
#' @return an age-period matrix od split population counts with the same number of 
#' columns as \code{popmat}, and single ages in rows.
#' 
#' @references 
#' Shryock, H. S., Siegel, J. S., & Larmon, E. A. (1973). 
#' The methods and materials of demography. US Bureau of the Census.

#' @export
#' 
#' @examples 
#' p5 <- structure(c(54170.08, 44774.6, 42141.587, 38463.515, 34405.607, 
#' 162369.816, 57424.3568738, 44475.4981681, 41751.7574114, 39628.4338929, 
#' 34756.9473002, 164194.0485702, 60272.2061248, 44780.1982856, 
#' 41803.6541424, 40229.0292664, 35154.7682192, 166275.9022992, 
#' 62726.896388, 45681.1355532, 42100.72506, 40473.8600572, 35598.545404, 
#' 168556.5331816, 64815.5458002, 47136.5341033, 42508.3026466, 
#' 40532.3096745, 36082.7490698, 170990.473735, 66579.122, 49070.407, 
#' 42953.604, 40534.586, 36596.844, 173545.633), .Dim = c(6L, 6L
#' ), .Dimnames = list(seq(0,25,5), 1950:1955))
#' head(p5) # this is the entire matrix
#' p1 <- grabill(p5)
#' head(p1); tail(p1)
#' colSums(p1) - colSums(p5) 
#' p1 - spragueSimple(p5)

grabill <- function(popmat){
	
	# get coefficient matrices for Sprague and Grabill
	scmg              <- grabillExpand(popmat)
	scm               <- spragueExpand(popmat)
	
	# split pop counts
	pops              <- scm %*% popmat
	popg              <- scmg %*% popmat
	
	# ---------------------------------------------
	# now we graft the two estimates in together,
	# preserving the middle part for grabill, and blending
	# aggressively into the young and closeout parts of Sprague
	# weights for grafting in grabill
	m                 <- nrow(pops)
	lr                <- m - 1
	fr                <- lr - 9
	
	# these weights do much better than linear weights.
	w10               <- exp(row(pops[1:10, ]) ) / exp(10.1)
	
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
	middle.sums       <- colSums(middle.part)
	# the difference to redistribute
	add.in            <- (middle.part %*% diag(1 / middle.sums)) %*% diag(redist)
	popg              <- popg + add.in
	# ---------------------------------------------
	# label dims and return
	rg                <- range(as.integer(rownames(popmat)))
	dimnames(popg)    <- list(rg[1]:rg[2], colnames(popmat))
	
	popg
}


#' split age groups using a monotonic spline
#' @description Take the cumulative sum of \code{Value} and then run a monotonic spline through it. The first 
#' differences split back single-age estimates of \code{Value}. Optionally keep the open age group untouched.
#' 
#' @param Value numeric vector of counts in age groups
#' @param Age5 integer vector of lower bound of age groups
#' @param keep.OAG logical (default \code{FALSE}). Would we like to re-impute the last 
#' element of \code{Value} as the open age group?
#' @return numeric vector of single age counts 
#' @export
#' @examples
#' Value <- structure(c(88623.0176512, 90841.8228447, 93438.8052066, 96324.9863902, 
#' 				99281.3083695, 102050.512557, 104351.4333985, 106554.7539415, 
#' 				109169.72444, 112188.3656672, 113582.4293037, 112613.9259593, 
#' 				108903.7255528, 102622.1778867, 95866.994701, 80874.3840953, 
#' 				60195.7933955, 37523.1859903, 17927.2905862, 5642.0540798, 1110.2324251
#' 		), .Names = c("0", "5", "10", "15", "20", "25", "30", "35", "40", 
#' 				"45", "50", "55", "60", "65", "70", "75", "80", "85", "90", "95", 
#' 				"100"))
#' 
#' splitMono(Value)
#' splitMono(Value, keep.OAG = TRUE)

splitMono <- function(Value, Age5 = seq(0, length(Value)*5-5, 5), keep.OAG = FALSE){
	AgePred    <- min(Age5):(max(Age5) + 1)
	y          <- c(0, cumsum(Value))
	x          <- c(0, Age5 + 5)
	y1         <- splinefun(y ~ x, method = "monoH.FC")(AgePred)
	single.out <- diff(y1)
	if (keep.OAG){
		single.out[length(single.out)] <- Value[length(Value)]
	}
	single.out
}

#' blend the Sprague upper boundary age estimates into monotonic spline estimates
#' 
#' @description A simple monotonic spline on the cumulative sum of population counts
#' may return more convincing single age count estimates than the Sprague splitting method.
#' This function blends the Sprague estimates starting at \code{pivotAge}.
#' 
#' @param popmat a numeric matrix of population counts in 5-year age groups, with integer-labeled 
#' margins (age in rows and year in columns).
#' @param pops optional numeric matrix of single age population counts derived from \code{popmat}.
#' @param pivotAge integer (default 90). Age at which to switch to spline-based estimates.
#' 
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
#'  popmat <- structure(c(54170.08, 44774.6, 42141.587, 38463.515, 34405.607, 
#' 30386.087, 26933.215, 23481.447, 20601.952, 16488.665, 14248.402, 
#' 9928.03, 8489.895, 4800.567, 3598.98, 2048.337, 940.945, 326.4, 
#' 80.117, 16.777, 0, 57424.3568738, 44475.4981681, 41751.7574114, 
#' 39628.4338929, 34756.9473002, 30605.4764074, 27182.9274184, 23792.1416427, 
#' 20723.6990619, 17056.1374644, 14058.6774471, 10585.1369045, 8102.9339783, 
#' 5306.0245132, 3366.8557838, 2040.3836514, 963.1478214, 314.5586137, 
#' 79.8410136, 15.502702, 0.6041464, 60272.2061248, 44780.1982856, 
#' 41803.6541424, 40229.0292664, 35154.7682192, 30977.7459584, 27456.1890384, 
#' 24096.7577352, 20872.8705064, 17545.6344304, 13989.6245496, 11145.817084, 
#' 7840.7544168, 5737.6492352, 3184.2782128, 2061.6884304, 960.8241584, 
#' 310.8441432, 79.5691376, 14.577848, 1.0774144, 62726.896388, 
#' 45681.1355532, 42100.72506, 40473.8600572, 35598.545404, 31438.8560456, 
#' 27758.1143968, 24396.368626, 21055.4431156, 17958.2723088, 14045.871538, 
#' 11588.5228268, 7731.3586756, 6059.5514704, 3085.9480552, 2082.5787944, 
#' 948.7749752, 312.0351852, 79.4346976, 13.97192, 1.4305504, 64815.5458002, 
#' 47136.5341033, 42508.3026466, 40532.3096745, 36082.7490698, 31939.7069818, 
#' 28091.9515352, 24693.2568243, 21274.4768803, 18299.2816212, 14223.4688239, 
#' 11906.0609901, 7785.1704959, 6254.5961836, 3090.4472686, 2083.8345794, 
#' 937.704655, 315.6524977, 79.5367336, 13.650846, 1.6768184), .Dim = c(21L, 
#' 5L), .Dimnames = list(c("0", "5", "10", "15", "20", "25", "30", 
#' "35", "40", "45", "50", "55", "60", "65", "70", "75", "80", "85", 
#' "90", "95", "100"), c("1950", "1951", "1952", "1953", "1954")))
#' closed.out <- spragueCloseout(popmat)
#' colSums(closed.out) - colSums(popmat)
#' spragueCloseout(popmat, pivotAge = 85)
#' # giving a different single-age split to close out this way:
#' popg <- grabill(popmat)
#' grabill.closed.out <- spragueCloseout(popmat, popg)
#' # totals not necessarily preserved if mixed w Grabill
#' # I wouldn't recommend a rescale of the total, since the 
#' # only part we mess with here is the old age section. Ergo,
#' # one may wish to instead rescale results colSums() of 
#' # popg at age pivotAge and higher.
#' colSums(grabill.closed.out) - colSums(popmat)

spragueCloseout <- function(popmat, pops, pivotAge = 90){
	if (missing(pops)){
		pops    <- spragueSimple(popmat)
	}
	# get the spline population split
	popmono <- apply(popmat, 2, splitMono, keep.OAG = TRUE)
	
	# some age pars
	Age5    <- as.integer(rownames(popmat))
	Age1    <- min(Age5):max(Age5)
	
	# some checks on pivotAge...
	if (!(max(Age1) - 10) >= pivotAge){
		pivotAge <- max(Age1) - 10
		if (pivotAge < 80){
			warning("pivotAge wasn't in rownames(popmat), moved it to 3rd 
from bottom row of popmat, but appears to be < 80
so returning spragueSimple() output as-is, no extra closeout performed.")
            return(pops)
		}
		warning("pivotAge moved to ", pivotAge, ", continued.")
	}
	# -----------------------------
	# now begin the closeout blend.
	p.i              <- which(Age1 == pivotAge)
	## substitute Sprague interpolation if > pchip for better belnding of the two series
	pop.c            <- popmono[p.i:(p.i + 4), ]
	ind              <- pops[p.i, ] > pop.c[1, ]
	pop.c[1, ind]    <- pops[p.i, ind] 

	## adjust back on initial pop 5x5 for age 90-94
	## proportional distribution
	pop.c[is.na(pop.c )] <- 0
	prop             <- prop.table(pop.c, margin = 2)
	pivot5           <- popmat[as.character(pivotAge), ]
    pop.c            <- prop %*% diag(pivot5)
	## append the remaining of the age groups (except last open age)
	## 95-99 onward
	m                <- nrow(pops)
	pop.c            <- rbind(pop.c, popmono[(p.i + 5):m, ])
	## append Sprague interpolation before age 90
	pop.c            <- rbind(pops[1:(p.i - 1), ], pop.c)
	
	## deal with negative values if applicable (but in principle should not be happening)
	pop.c[pop.c < 0] <- 0
	
	# label and return
	dimnames(pop.c) <- list(Age1, colnames(popmat))
	
	pop.c
}
