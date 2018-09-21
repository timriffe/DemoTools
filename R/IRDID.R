
# Author: tim
###############################################################################

#' Index of relative difference.
#' 
#' @description Calculate the relative percent difference between two population structures. A returned 
#' value of zero means that the two population have identical structure. 
#' 
#' @param pop1 numeric. Vector of counts from population 1.
#' @param pop2 numeric. Vector of counts from population 2.
#' 
#' @details Input populations are assumed to be ordered in the same way prior to calling the function. 
#' It is only checked  that the vectors are of the same length. The input arguments could indeed be populations structured on multiple
#' variables (more than just age), as long as they are ordered in the same way. It is advised to lower 
#' the open age group for this method because each age has the same weight. Ages where one population 
#' has a zero count and the other does not are thrown out.
#' 
#' @return The value of the index ranging from 0 to infinity.
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
	n     <- length(pop1)
	
	r1    <- rescale.vector(pop1)
	r2    <- rescale.vector(pop2)
	
	ratio <- r2 / r1 
	# we don't like zeros in denominators
	ratio[is.infinite(ratio)] <- NA
	sum(abs((ratio * 100) - 100), na.rm = TRUE) / (2 * n)
	
}

#' Index of dissimilarity
#' 
#' @description Calculate the index of dissimilarity between two population structures. A returned 
#' value of zero means that the two population have identical structure. A value of 100 means
#' that the populations have no overlap at all (not likely for populations structured only by age). 
#' This is a simple measure of distribution overlap on the absolute scale.
#' 
#' @details Input populations are assumed to be ordered in the same way prior to calling the function. 
#' It is only checked  that the vectors are of the same length. The input arguments could indeed be populations structured on multiple
#' variables (more than just age), as long as they are ordered in the same way.
#' 
#' @param pop1 numeric. Vector of counts from population 1.
#' @param pop2 numeric. Vector of counts from population 2.
#' 
#' @references 
#' \insertRef{siegel2004methods}{DemoTools}
#' @return The value of the index ranging from 0 to infinity..
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

