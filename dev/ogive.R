
# Author: tim
###############################################################################

# notes from patrick:
# should be a generic function with loess as first option under the hood,
# choice between user-specified smoothing paramater and gcv crit
# 2nd option would be 2nd degree polynomial (no auto detect, but in general
# intended when there are few points, like fert in 5-yr age groups), but don't sweat too much
# about endpoints as there are purpose-designed methods for particular 
# phenomena.


# example data taken from ZELNIK.R
Value <- c(941307,1041335,1237034,1411359,1383853,1541942,1321576,1285877,1563448,886705,
       1623998,562924,1485173,543216,771219,903496,686431,370007,942999,250820,
       1023667,200131,688640,222011,281738,1239965,288363,263326,483143,78635,
       1349886,68438,415127,101596,100758,1392434,178633,126351,286520,50836,
       1331036,48995,251153,58393,54995,1033812,68792,72766,175943,28254,
       1038747,32894,136179,37667,38230,596049,52602,36493,74106,16759,
       790643,20596,70109,18044,19891,357491,15253,17489,31057,8481,
       429816,7951,35583,8612,6589,454645)
Age  <- 0:75

# ---------
# less obvious arguments:
# OAG: logical is the last value an open age group that shoud be preserved?
# ... optional arguments to pass to \code{loess()}. notably \code{span} might be interesting, as it controls smoothness
ogive_loess <- function(Value, Age, OAG = TRUE, ...){
	scale      <- sum(Value, na.rm = TRUE)
	N          <- length(Value)
	stopifnot(N == length(Age))
	
	age        <- Age
	
	# separate Open age group if desired.
	if (OAG){
		OA     <- Value[N]
		Value  <- Value[-N]
		age    <- age[-N]
		scale  <- scale - OA
	}
	
	fit        <- loess(Value ~ age, ...)
	out        <- fit$fitted
	
	# enforce sums
	out        <- rescale.vector(out,  scale) 
	
	# tack open age group back on
	if (OAG){
		out    <- c(out, OA)
	}
	
	# name vector
	names(out) <- Age
	
	out
}


# ---------
# less obvious arguments:
# degree integer degree of polynomial, default 2
# trans character. tranformation to Value prior to fitting? \code{"log"} is the
# only implemented possibility at this time. Leave empty otherwise.
# OAG: logical is the last value an open age group that shoud be preserved?
# ... other optional arguments to pass to \code{lm()}
ogive_poly <- function(Value, Age, degree = 2, trans, OAG = TRUE, ...){
	scale     <- sum(Value, na.rm = TRUE)
	N         <- length(Value)
	stopifnot(N == length(Age))
	
	# rename because we modify
	age       <- Age
	value     <- Value
	
	# default regression weights
	w         <- rep(1, length(value))
	# separate Open age group if desired.
	if (OAG){
		OA    <- value[N]
		value <- value[-N]
		age   <- age[-N]
		scale <- scale - OA
		w     <- w[-N]
	}
	
	# -------------------------
	# log transform to constrain positive, if desired
	# indicator
	lTF       <- FALSE
	if (!missing(trans)){
		if (trans == "log"){
			lTF   <- TRUE
			value <- log(value)
			w[is.infinite(value)] <- 0
		} else {
			warning("The only transformation implemented now is log(Value)\nSpecify trans = 'log' for this.\nfunction continued with no transformation.")
		}
		
	}
	
	# -------------------------
	# build up polynomial expression
	polys     <- 2:degree
	polys2    <- paste0("age^", polys)
	polys3    <- paste0(" + I(",polys2,")")
	# final formula
	expr      <- paste0("Value. ~ age", paste(polys3, collapse = ""))
	
	# stick elements into ad hoc data.frame
	dataf     <- data.frame(Value. = value, age = age, w = w)
	# fit linear model
	fit       <- lm(expr, weight = w, data = dataf, ...)
	# get predicted values
	out       <- predict(fit)
	
	# -------------------------
	# transform back, if necessary
	if (lTF){
			out <- exp(out)
	} 		
	
	# -------------------------
    # closing steps
    # enforce sums
    out        <- rescale.vector(out,  scale) 
	
	# tack OAG back on, if necessary
	if (OAG){
		out    <- c(out, OA)
	}
	
	# name vector
	names(out) <- Age
	
	out
}


# same args as above:

ogive <- function(Value, Age, method = "loess", OAG = TRUE, degree = 2, trans, ...){
	
	if (method == "loess"){
		out <- ogive_loess(Value = Value, Age = Age, OAG = OAG, ...)
	}
	if (method == "poly"){
		out <- ogive_poly(Value = Value, Age = Age, degree = degree, trans = trans, OAG = OAG, ...)
	}
	out
}


# examples
# \dontrun{
#plot(Age, Value, pch = 16, col = gray(.5))
#lines(Age, ogive_poly(Value, Age, degree = 2))
#lines(Age, ogive_poly(Value, Age, degree = 3))
#lines(Age, ogive_poly(Value, Age, degree = 2, OAG = FALSE))
#lines(Age, ogive_poly(Value, Age, degree = 3, OAG = FALSE))
#
#lines(Age, ogive_poly(Value, Age, degree = 2, trans = "log"))
#lines(Age, ogive_poly(Value, Age, degree = 3, trans = "log"))
#lines(Age, ogive_poly(Value, Age, degree = 2, OAG = FALSE, trans = "log"))
#lines(Age, ogive_poly(Value, Age, degree = 3, OAG = FALSE, trans = "log"))

# different only toward the upper age boundary due to OAG
#plot(Age, Value, pch = 16, col = gray(.5))
#lines(Age, ogive_loess(Value, Age, OAG = FALSE), col = "blue",lwd = 2)
#lines(Age, ogive_loess(Value, Age, OAG = TRUE), col = "red",lwd = 2)
## you can change the span argument for loess to modify smoothness.
#lines(Age, ogive_loess(Value, Age, OAG = FALSE, span = .5), col = "blue",lwd = 2, lty = 2)
#lines(Age, ogive_loess(Value, Age, OAG = FALSE, span = 2), col = "blue",lwd = 2, lty = 2)
#
#}