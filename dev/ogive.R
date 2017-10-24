
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


ogive_loess <- function(Value, Age, scale = sum(Value), OAG = TRUE, ...){
	
	N   <- length(Value)
	stopifnot(N == length(Age))
	
	age <- Age
	
	# separate Open age group if desired.
	if (OAG){
		OA    <- Value[N]
		Value <- Value[-N]
		age   <- age[-N]
	}
	
	fit <- loess(Value ~ age, ...)
	out <- fit$fitted
	
	if (OAG){
		out <- c(out, OA)
	}
	
	names(out) <- Age
	
	# enforce sums
	out <- rescale.vector(out) * scale
	
	out
}

# TODO log transform not working...
ogive_poly <- function(Value, Age, degree = 2, trans, scale = sum(Value), OAG = TRUE, ...){

	N         <- length(Value)
	stopifnot(N == length(Age))
	
	age       <- Age
	
	# separate Open age group if desired.
	if (OAG){
		OA    <- Value[N]
		Value <- Value[-N]
		age   <- age[-N]
	}
	
	
	# log transform to constrain positive, if desired
	w     <- rep(1, length(Value))
	if (!missing(trans)){
		if (trans == "log"){
			lTF <- TRUE
			Value <- log(Value)
		} else {
			warning("The only transformation implemented now is log(Value)\nSpecify trans = 'log' for this.\nfunction continued with no transformation.")
		}
		
	}
	
	polys     <- 2:degree
	polys2    <- paste0("age^", polys)
	polys3    <- paste0(" + I(",polys2,")")
	expr      <- paste0("Value ~ age", paste(polys3, collapse = ""))
	data      <- data.frame(Value = Value, age = age, w = w)
	fit       <- lm(expr, weight = w, data = data)
	out       <- predict(fit)
	
	# transform back, if necessary
	if (lTF){
			out <- exp(out)
	} 		
	
	
	# enforce sums
	out        <- rescale.vector(out) * scale
	
	# tack OAG back on, if necessary
	if (OAG){
		out    <- c(out, OA)
	}
	names(out) <- Age
	
	out
}

plot(Age, Value, pch = 16, col = gray(.5))
lines(Age, ogive_poly(Value, Age, degree = 2))
lines(Age, ogive_poly(Value, Age, degree = 3))
lines(Age, ogive_poly(Value, Age, degree = 2, OAG = FALSE))
lines(Age, ogive_poly(Value, Age, degree = 3, OAG = FALSE))

lines(Age, ogive_poly(Value, Age, degree = 2, trans = "log"))
lines(Age, ogive_poly(Value, Age, degree = 3))
lines(Age, ogive_poly(Value, Age, degree = 2, OAG = FALSE))
lines(Age, ogive_poly(Value, Age, degree = 3, OAG = FALSE))