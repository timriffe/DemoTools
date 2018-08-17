
# Author: tim
###############################################################################
# TODO: 1) splitMono() should be flexible wrt AgeInt, perhaps taking an AgeInt arg.
# follow same detection protocol as in splitUniform to produce intermediate pieces 
# needed. Be sure to extend over final interval if OAG = FALSE.

# Counts -> Value
splitUniform <- function(Value, AgeInt, Age, OAG = TRUE, OAvalue = 1){
	if (missing(Age) & missing(AgeInt)){
		Age <- names2age(Value)
	}
	if (missing(AgeInt)){
		# give 1 to final interval to preserve
		AgeInt <- age2int(Age, OAG = OAG, OAvalue = OAvalue)
	}
	if (missing(Age)){
		Age <- cumsum(AgeInt) - AgeInt
	}
	# discount for single
	out <- rep(Value / AgeInt, times = AgeInt)
	names(out) <- min(Age):(min(Age) + length(out) - 1)
	out
}

# Prob: splitMono being used as Mono_smth() would have been intended.
# make splitting and smoothing separate so avoid confusion.

# we want splitMono() to emulate this behavior
# monoCloseout() depends on this, and splitOscillate() depends on monoCloseout().
# be sure to be able to get present dimensions of splitMono() outpout as option.
splitMono <- function(Value, AgeInt, Age, OAG = FALSE){
	if (missing(Age) & missing(AgeInt)){
		Age    <- names2age(Value)
	}
	if (missing(AgeInt)){
		# give 1 to final interval to preserve
		AgeInt <- age2int(Age, OAG = OAG, OAvalue = 1)
	}
	if (missing(Age)){
		Age    <- int2age(AgeInt)
	}

    # if age is single return as-is
	if (is_single(Age)){
		return(Value)
	}
	
	# if last age Open, we preserve it
	if (OAG){
		N       <- length(Value)
		OAvalue <- Value[N]
		Value   <- Value[-N]
		Age     <- Age[-N]
		AgeInt  <- AgeInt[-N]
	}
	# if the final age is Open, then we should remove it and then
	# stick it back on
	
	AgePred    <- c(min(Age), cumsum(AgeInt))
	y          <- c(0, cumsum(Value))
	AgeS       <- min(Age):sum(AgeInt) 
	y1         <- splinefun(y ~ AgePred, method = "monoH.FC")(AgeS)
	out        <- diff(y1)
	names(out) <- AgeS[-length(AgeS)]
 
	# The open age group is maintained as-is.
	if (OAG){
		out              <- c(out, OAvalue)
		names(out)       <- AgeS
	}
	
	out
}



splitMono2(Value, OAG = FALSE)      # OAG = TRUE implies identical behavior,
splitUniform(Value, OAG = FALSE)   # OAG = FALSE splitUniform extends final age group through,

plot(0:104, splitUniform(Value, OAG = FALSE),type='s')
lines(0:104, splitMono2(Value, OAG = FALSE))

plot(0:100, splitUniform(Value, OAG = TRUE),type='s')
lines(0:100, splitMono2(Value, OAG = TRUE))

# if OAG = FALSE, identical shape

 Value <- structure(c(88623, 90842, 93439, 96325, 99281, 102051, 104351, 
				 106555, 109170, 112188, 113582, 112614, 108904, 102622, 95867, 
				 80874, 60196, 37523, 17927, 5642, 1110), .Names = c("0", "5", 
				 "10", "15", "20", "25", "30", "35", "40", "45", "50", "55", "60", 
				 "65", "70", "75", "80", "85", "90", "95", "100"))
 
 
 
?rescaleAgeGroups

set.seed(3)
AgeIntRandom <- sample(1:5, size = 15,replace = TRUE)
AgeInt5      <- rep(5, 9)
original     <- runif(45, min = 0, max = 100)
pop1         <- groupAges(original, 0:45, AgeN = int2ageN(AgeIntRandom, FALSE))
pop2         <- groupAges(original, 0:45, AgeN = int2ageN(AgeInt5, FALSE))
# inflate (in this case) pop2
perturb      <- runif(length(pop2), min = 1.05, max = 1.2)

pop2         <- pop2 * perturb

# a recursively constrained solution
(sMonRF <- rescaleAgeGroups(Value1 = pop1, 
					AgeInt1 = AgeIntRandom,
					Value2 = pop2, 
					AgeInt2 = AgeInt5, 
					splitfun = splitMono, 
					recursive = FALSE))
(sMonRT <- rescaleAgeGroups(Value1 = pop1, 
					AgeInt1 = AgeIntRandom,
					Value2 = pop2, 
					AgeInt2 = AgeInt5, 
					splitfun = splitMono, 
					recursive = TRUE))

sMonRF - sMonRT




popmat <- structure(c(54170, 44775, 42142, 38464, 34406, 30386, 26933, 
				23481, 20602, 16489, 14248, 9928, 8490, 4801, 3599, 2048, 941, 
				326, 80, 17, 0, 57424, 44475, 41752, 39628, 34757, 30605, 27183, 
				23792, 20724, 17056, 14059, 10585, 8103, 5306, 3367, 2040, 963, 
				315, 80, 16, 1, 60272, 44780, 41804, 40229, 35155, 30978, 27456, 
				24097, 20873, 17546, 13990, 11146, 7841, 5738, 3184, 2062, 961, 
				311, 80, 15, 1, 62727, 45681, 42101, 40474, 35599, 31439, 27758, 
				24396, 21055, 17958, 14046, 11589, 7731, 6060, 3086, 2083, 949, 
				312, 79, 14, 1, 64816, 47137, 42508, 40532, 36083, 31940, 28092, 
				24693, 21274, 18299, 14223, 11906, 7785, 6255, 3090, 2084, 938, 
				316, 80, 14, 2), .Dim = c(21L, 5L), .Dimnames = list(c("0", "5", 
						"10", "15", "20", "25", "30", "35", "40", "45", "50", "55", "60", 
						"65", "70", "75", "80", "85", "90", "95", "100"), c("1950", "1951", 
						"1952", "1953", "1954")))


set.seed(3)
AgeIntRandom <- sample(1:5, size = 15,replace = TRUE)
AgeInt5      <- rep(5, 9)
original     <- runif(45, min = 0, max = 100)
pop1         <- groupAges(original, 0:45, AgeN = int2ageN(AgeIntRandom, FALSE))
pop2         <- groupAges(original, 0:45, AgeN = int2ageN(AgeInt5, FALSE))
# inflate (in this case) pop2
perturb      <- runif(length(pop2), min = 1.05, max = 1.2)

pop2         <- pop2 * perturb

# a recursively constrained solution
(pop1resc <- rescaleAgeGroups(Value1 = pop1, 
					AgeInt1 = AgeIntRandom,
					Value2 = pop2, 
					AgeInt2 = AgeInt5, 
					splitfun = splitUniform, 
					recursive = TRUE))
# a single pass adjustment (no recursion)
(pop1resc1 <- rescaleAgeGroups(Value1 = pop1, 
					AgeInt1 = AgeIntRandom,
					Value2 = pop2, 
					AgeInt2 = AgeInt5, 
					splitfun = splitUniform, 
					recursive = FALSE))
(pop1rescm <- rescaleAgeGroups(Value1 = pop1, 
					AgeInt1 = AgeIntRandom,
					Value2 = pop2, 
					AgeInt2 = AgeInt5, 
					splitfun = splitMono, 
					recursive = TRUE))
# a single pass adjustment (no recursion)
(pop1resc1m <- rescaleAgeGroups(Value1 = pop1, 
					AgeInt1 = AgeIntRandom,
					Value2 = pop2, 
					AgeInt2 = AgeInt5, 
					splitfun = splitMono, 
					recursive = FALSE))



## Not run:


# show before / after
plot(NULL,xlim=c(0,45),ylim=c(0,2),main = "Different (but integer) intervals",
		xlab = "Age", ylab = "", axes = FALSE)
x1 <- c(0,cumsum(AgeIntRandom))
rect(x1[-length(x1)],1,x1[-1],2,col = gray(.8), border = "white")
x2 <- c(0,cumsum(AgeInt5))
rect(x2[-length(x2)],0,x2[-1],1,col = "palegreen1", border = "white")
text(23,1.5,"Original (arbitrary grouping)",font = 2, cex=1.5)
text(23,.5,"Standard to rescale to (arbitrary grouping)",font = 2, cex=1.5)
axis(1)

# adjustment factors:
plot(int2age(AgeInt5), perturb, ylim = c(0, 2))
points(int2age(AgeIntRandom), pop1resc / pop1, pch = 16)
# non-recursive is less disruptive for uniform
points(int2age(AgeIntRandom), pop1resc1 / pop1, pch = 16, col = "blue")

# show before / after under uniform (in pop1) assumption.
plot(NULL, xlim = c(0, 45), ylim = c(0, 150), main = "Uniform constraint")
lines(0:44, splitUniform(pop1, AgeInt = AgeIntRandom, OAG = FALSE), col = "red")
lines(0:44, splitUniform(pop2, AgeInt = AgeInt5, OAG = FALSE), col = "blue")
lines(0:44, splitUniform(pop1resc, AgeInt = AgeIntRandom, OAG = FALSE), 
		col = "orange", lty = 2, lwd = 2)
lines(0:44, splitUniform(pop1resc1, AgeInt = AgeIntRandom, OAG = FALSE), 
		col = "magenta", lty = 2, lwd = 2)
legend("topright",
		lty = c(1, 1, 2, 2),
		col = c("red", "blue", "orange", "magenta"),
		lwd = c(1, 1, 2, 2),
		legend = c("Original N1", "Prior N2",
				"Rescaled N1 recursive", "Rescaled N1 1 pass"))
## End(Not run)

plot(NULL, xlim = c(0, 45), ylim = c(0, 150), main = "Uniform constraint")
lines(0:44, splitUniform(pop1, AgeInt = AgeIntRandom, OAG = FALSE), col = "red")
lines(0:44, splitUniform(pop2, AgeInt = AgeInt5, OAG = FALSE), col = "blue")

lines(0:44, splitUniform(pop1resc, AgeInt = AgeIntRandom, OAG = FALSE), 
		col = "orange", lty = 2, lwd = 2)
lines(0:44, splitUniform(pop1resc1, AgeInt = AgeIntRandom, OAG = FALSE), 
		col = "magenta", lty = 2, lwd = 2)


#legend("topright",
#		lty = c(1, 1, 2, 2),
#		col = c("red", "blue", "orange", "magenta"),
#		lwd = c(1, 1, 2, 2),
#		legend = c("Original N1", "Prior N2",
#				"Rescaled N1 recursive", "Rescaled N1 1 pass"))


plot(NULL, xlim = c(0, 45), ylim = c(0, 150), main = "Uniform constraint")
lines(0:44, splitMono(pop1, AgeInt = AgeIntRandom, OAG = FALSE), col = "red")
lines(0:44, splitMono(pop2, AgeInt = AgeInt5, OAG = FALSE), col = "blue")

lines(0:44, splitMono(pop1rescm, AgeInt = AgeIntRandom, OAG = FALSE), 
		col = "orange", lty = 2, lwd = 2)
lines(0:44, splitMono(pop1resc1m, AgeInt = AgeIntRandom, OAG = FALSE), 
		col = "magenta", lty = 2, lwd = 2)