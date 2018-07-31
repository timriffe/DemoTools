
# Author: tim
###############################################################################

# PAS BASEPOP spreadsheet.

# census date
ReferenceDate <- 1986.21 


date1 <- 1977.31
date2 <- 1986.50 
Lx1 <- c(87732, 304435, 361064)
Lx2 <- c(88451, 310605, 370362)

# from spreadsheet E64-E66
result1test <- c(88389, 310075, 369563)
result2test <- c(88233, 308732, 367539)
result3test <- c(87842, 305375, 362480)
result1 <- interpolatePop(
		Pop1 = Lx1, 
		Pop2 = Lx2, 
		Date1 = date1, 
		Date2 = date2, 
		DesiredDate = ReferenceDate - .5, 
		method = "linear")
result2 <- interpolatePop(
		Pop1 = Lx1, 
		Pop2 = Lx2, 
		Date1 = date1, 
		Date2 = date2, 
		DesiredDate = ReferenceDate - 2.5, 
		method = "linear")
result3 <- interpolatePop(
		Pop1 = Lx1, 
		Pop2 = Lx2, 
		Date1 = date1, 
		Date2 = date2, 
		DesiredDate = ReferenceDate - 7.5, 
		method = "linear")
stopifnot(all((round(result1) - result1test) == 0))
stopifnot(all((round(result2) - result2test) == 0))
stopifnot(all((round(result2) - result2test) == 0))

# but can it extrapolate if necessary?
# this will work just the same for asfr. Would like to allow for arbitrary number of 
# input Lx vectors. Also ideally in varying shape and size. can do 2d linear approx:
# first over age to get matching ages, then over time to get desired ref points.




