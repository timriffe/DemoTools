
# Author: tim
###############################################################################

#' the basic Sprague age-splitting method
#' 
#' @description This method is based on the first stage of the Sprague R 
#' script prepared by Thomas Buettner and Patrick Gerland, utself based on the description
#' in Siegel and Swanson, 2004, p. 727. This is based on a set of prior fixed coefficients, not a spline.
#' 
#' @param popmat a numeric matrix of population counts in 5-year age groups, with integer-labeled margins (age in rows and year in columns).
#' 
#' @details Ages should refer to lower age bounds, ending in the open age group (not a closed terminal age). 
#' It is important that \code{popmat} is a \code{matrix} and not a \code{data.frame}, and the dimension labelling
#' is also necessary.
#' 
#' @return an age-period matrix od split population counts with the same number of 
#' columns as \code{popmat}, and single ages in rows.
#' 
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
#' colSums(p1) - colSums(p5) # TODO fix bug. colSums(scm) not all 1.
spragueSimple <- function(popmat){
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
	
	## Sprague's Split
	bc <- c(
			0.3616, -0.2768,  0.1488, -0.0336,  0.0000, 
			0.2640, -0.0960,  0.0400, -0.0080,  0.0000, 
			0.1840,  0.0400, -0.0320,  0.0080,  0.0000, 
			0.1200,  0.1360, -0.0720,  0.0160,  0.0000, 
			0.0704,  0.1968, -0.0848,  0.0176,  0.0000, 
			0.0336,  0.2272, -0.0752,  0.0144,  0.0000, 
			0.0080,  0.2320, -0.0480,  0.0080,  0.0000, 
			-0.0080,  0.2160, -0.0080,  0.0000,  0.0000, 
			-0.0160,  0.1840,  0.0400, -0.0080,  0.0000, 
			-0.0176,  0.1408,  0.0912, -0.0144,  0.0000, 
			-0.0128,  0.0848,  0.1504, -0.0240,  0.0016, 
			-0.0016,  0.0144,  0.2224, -0.0416,  0.0064, 
			0.0064, -0.0336,  0.2544, -0.0336,  0.0064, 
			0.0064, -0.0416,  0.2224,  0.0144, -0.0016, 
			0.0016, -0.0240,  0.1504,  0.0848, -0.0128, 
			0.0000, -0.0144,  0.0912,  0.1408, -0.0176, 
			0.0000, -0.0080,  0.0400,  0.1840, -0.0160, 
			0.0000,  0.0000, -0.0080,  0.2160, -0.0080, 
			0.0000,  0.0080, -0.0480,  0.2320,  0.0080, 
			0.0000,  0.0144, -0.0752,  0.2272,  0.0336, 
			0.0000,  0.0176, -0.0848,  0.1968,  0.0704, 
			0.0000,  0.0160, -0.0720,  0.1360,  0.1200, 
			0.0000,  0.0080, -0.0320,  0.0400,  0.1840, 
			0.0000, -0.0080,  0.0400, -0.0960,  0.2640, 
			0.0000, -0.0336,  0.1488, -0.2768,  0.3616
	)
	## standard format for Sprague coefficients
	sm             <- matrix(bc, 25, 5, byrow = TRUE)
	
	## create a Sprague coefficient matrix for 5-year age groups
	scm            <- array(0, dim = c(m1 + 1, m))
	
	## insert upper left block
	scm[1:10, 1:5] <- sm[1:10, ]
	
	# determine positions of middle blocks
	rowpos         <- matrix(10:((MP*5) + 9), ncol = 5, byrow = TRUE)
	colpos         <- row(rowpos) + col(rowpos) - 1
	for (i in (1:MP)) {
		# calculate the slices and add middle panels accordingly
		scm[rowpos[i,], colpos[i, ]] <- sm[11:15,]
	}
	
	## insert last two panels
	fr                <- (m - 3) * 5 + 1
	lr                <- (m - 1) * 5
	fc                <- MP 
	lc                <- MP + 4 
	scm[fr:lr,fc:lc]  <- sm[16:25,]
	
	# last open ended age group
	scm[m1 + 1, m]    <- 1
	
	pop1              <- scm %*% popmat
	dimnames(pop1) <- list(Age1, yrs)
	pop1
}


