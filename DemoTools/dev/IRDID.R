
# Author: tim
###############################################################################

#' index of relative difference
#' 
#' @description Calculate the relative percent difference between two population structures. A returned 
#' value of zero means that the two population have identical structure. A value of 100 means
#' that the populations have no overlap at all (not likely for populations structured only by age).
#' 
#' @details We assume that the populations are ordered in the same way and 
#' that said ordering happens before this function is called. We only check here that the vectors
#' are of the same length. The input arguments could indeed be popualtions structured on multiple
#' variables (more than just age), as long as they are ordered in the same way. It is advised to lower 
#' the open age group for this method because each age has the same weight.

#' @param pop1 numeric vector of counts from population 1
#' @param pop2 numeric vector of counts from population 2
#' 
#' @return an index value ranging from 0 to infinity.
#' @export
#' @examples
#' pop1 <- c(7.38,14.16,14.79,17.36,15.11,10.14,8.50,7.28,5.28)
#' pop2 <- c(6.48,12.27,15.25,15.10,14.66,10.80,8.95,9.28,7.21)
#' IRD(pop1, pop2)     # 6.7 reproduces table 7.20 of Siegel & Swanson (2004)
#' IRD(pop1,pop1)      # identical pops  = 0
#' IRD(pop1, pop2 * 2) # only structure matters
#' pop3      <- pop1
#' pop3[1:7] <- 0
#' IRD(pop1, pop3)     # theoretical max > 100 
IRD <- function(pop1, pop2){
	
	# for now just make sure they're the same lengths. We won't do
	# age matching here.
	stopifnot(length(pop1) == length(pop2))
	n  <- length(pop1)
	
	r1 <- rescale.vector(pop1)
	r2 <- rescale.vector(pop2)
	
	sum(abs((r2 / r1 * 100) - 100)) / (2 * n)
	
}

#' index of dissimilarity
#' 
#' @description Calculate the index of dissimilarity between two population structures. A returned 
#' value of zero means that the two population have identical structure. A value of 100 means
#' that the populations have no overlap at all (not likely for populations structured only by age). 
#' This is a simple measure of distribution overlap on the absolute scale.
#' 
#' @details We assume that the populations are ordered in the same way and 
#' that said ordering happens before this function is called. We only check here that the vectors
#' are of the same length. The input arguments could indeed be popualtions structured on multiple
#' variables (more than just age), as long as they are ordered in the same way.

#' @param pop1 numeric vector of counts from population 1
#' @param pop2 numeric vector of counts from population 2
#' 
#' @return an index value ranging from 0 to 100.
#' @export
#' @examples
#' pop1 <- c(7.38,14.16,14.79,17.36,15.11,10.14,8.50,7.28,5.28)
#' pop2 <- c(6.48,12.27,15.25,15.10,14.66,10.80,8.95,9.28,7.21)
#' ID(pop1, pop2)     # 5.5 reproduces table 7.20 of Siegel & Swanson (2004)
#' ID(pop1, pop1)     # identical = 0
#' ID(pop1, pop2 * 2) # scale invariant
#' pop3 <- pop4         <- pop1
#' pop3[1:5]            <- 0
#' pop4[6:length(pop4)] <- 0
#' ID(pop3, pop4)     # no overlap = 100.
ID <- function(pop1, pop2){
	# for now just make sure they're the same lengths. We won't do
	# age matching here.
	stopifnot(length(pop1) == length(pop2))
	n  <- length(pop1)
	
	r1 <- rescale.vector(pop1, 100)
	r2 <- rescale.vector(pop2, 100)
	
	sum(abs(r1 - r2)) / 2
}

