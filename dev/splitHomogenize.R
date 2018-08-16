
# Author: tim
###############################################################################


splitUniform <- function(Counts, AgeInt, Age, OAG = TRUE, OAvalue = 1){
	if (missing(AgeInt)){
		# give 1 to final interval to preserve
		AgeInt <- age2int(Age, OAG = OAG, OAvalue = OAvalue)
	}
	if (missing(Age)){
		Age <- cumsum(AgeInt) - AgeInt
	}
	# discount for single
	out <- rep(Counts / AgeInt, times = AgeInt)
	names(out) <- min(Age):(min(Age) + length(out) - 1)
	out
}


# we want splitMono() to emulate this behavior

splitMono <- function(Value, Age, OAG = FALSE){
	if (missing(Age)){
		Age               <- as.integer(names(Value))
	}
	
	# this is innocuous if ages are already grouped
	Age5       <- calcAgeN(Age, N = 5)
	pop5       <- groupAges(Value, Age = Age, AgeN = Age5, shiftdown = 0)
	
	AgePred    <- min(Age):(max(Age) + 1)
	y          <- c(0, cumsum(pop5))
	x          <- c(0, sort(unique(Age5)) + 5)
	y1         <- splinefun(y ~ x, method = "monoH.FC")(AgePred)
	single.out <- diff(y1)
	if (OAG){
		single.out[length(single.out)] <- pop5[length(pop5)]
	}
	
	names(single.out) <- min(Age):max(Age)
	single.out
}

