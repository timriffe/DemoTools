
# Author: tim
###############################################################################
# TODO: 1) splitMono() should be flexible wrt AgeInt, perhaps taking an AgeInt arg.
# follow same detection protocol as in splitUniform to produce intermediate pieces 
# needed. Be sure to extend over final interval if OAG = FALSE.

splitUniform <- function(Counts, AgeInt, Age, OAG = TRUE, OAvalue = 1){
	if (missing(Age) & missing(AgeInt)){
		Age <- names2age(Counts)
	}
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
# monoCloseout() depends on this, and splitOscillate() depends on monoCloseout().
# be sure to be able to get present dimensions of splitMono() outpout as option.
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
	
	# that is, the open age group is maintained as-is.
	if (OAG){
		single.out[length(single.out)] <- pop5[length(pop5)]
	}
	
	names(single.out) <- min(Age):max(Age)
	single.out
}
splitMono(Value, OAG = FALSE)      # OAG = TRUE implies identical behavior,
splitUniform(Value, OAG = FALSE)   # OAG = FALSE splitUniform extends final age group through,

# if OAG = FALSE, identical shape

 Value <- structure(c(88623, 90842, 93439, 96325, 99281, 102051, 104351, 
				 106555, 109170, 112188, 113582, 112614, 108904, 102622, 95867, 
				 80874, 60196, 37523, 17927, 5642, 1110), .Names = c("0", "5", 
				 "10", "15", "20", "25", "30", "35", "40", "45", "50", "55", "60", 
				 "65", "70", "75", "80", "85", "90", "95", "100"))