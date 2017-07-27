
# Author: tim
###############################################################################


# want some data w age heaping here:
Value <- c(80626,95823,104315,115813,100796,105086,97266,116328,
		75984,89982,95525,56973,78767,65672,53438,85014,
		47600,64363,42195,42262,73221,30080,34391,29072,
		20531,66171,24029,44227,24128,23599,82088,16454,
		22628,17108,12531,57325,17220,28425,16206,17532,
		65976,11593,15828,13541,8133,44696,11165,18543,
		12614,12041,55798,9324,10772,10453,6773,28358,
		9916,13348,8039,7583,42470,5288,5317,6582,
		3361,17949,3650,5873,3279,3336,27368,1965,
		2495,2319,1335,12022,1401,1668,1360,1185,
		9167,424,568,462,282,6206,343,409,333,291,4137,133,169,157,89,2068,68,81,66,57)
Age <- 0:99
#barplot(Value, main = "yup, these have heaping!")

Age0_5 <- calcAgeN(Age)
Age1_5 <- calcAgeN(Age, shiftdown = 1)
Age2_5 <- calcAgeN(Age, shiftdown = 2)
Age3_5 <- calcAgeN(Age, shiftdown = 3)
Age4_5 <- calcAgeN(Age, shiftdown = 4)

# select which ages to keep:
p1x1   <- matrix(nrow = length(Value), ncol = 5)
for (i in 0:4){
	Age.i.5  <- calcAgeN(Age, shiftdown = i)
	keep.i   <- rep(rle(Age.i.5)$leng, rle(Age.i.5)$leng) == 5
	Age.i    <- Age[keep.i]
	Val.i    <- Value[keep.i]
	Val.i.5  <- groupAges(Val.i, Age.i, AgeN = Age.i.5)
	# make fake open age
	Val.i.5  <- c(Val.i.5, pi)
	names(Val.i.5) <- c(unique(Age.i.5), max(Age.i.5) + 5)
	pop.est  <- spragueSimple(Val.i.5)
	p1x1[i + 1, keep.i] <- pop.est[as.character(Age.i), ]
}
keep.i <- rep(rle(Age1_5)$leng, rle(Age1_5)$leng) == 5



