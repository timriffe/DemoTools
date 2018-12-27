# Author: Tim Riffe
# Dec 27, 2017.. hahaha and Dec 27 2018 too!
# OPAG fresh from scratch, using mostly verbal description of how it should work.
# may or may not have a standard to test against, that will come later if so.
# ----------------------------------
# standard pop is a vector here. Could come from anywhere...
# make utilities to produce standards.
# assume single ages. can always graduate or group as required.
#' redistripute an open age group count over higher ages proportional to an arbitrary standard
#' @description This method could be useful whenever a reasonable standard is available. At present the standard must be supplied by the user.
#' @details In this implementation both the original population counts and the standard must be in single ages.
#' @param Pop numeric vector of population counts
#' @param Age integer vector of single age lower bounds
#' @param OAnow integer. The lower age bound above which counts will be redistributed
#' @param StPop numeric vector of standard population counts
#' @param StAge integer vector of single age lower bounds for the standard population
#' @param OAnew integer. The desired new open age, must be no higher than \code{max(StAge)}.
#' @export
#' @references 
#' \insertRef{PAS}{DemoTools}
#' @examples
#'  Pop <- c(38129,38382,38824,39275,39500,37304,35152,
#'  34061,33911,32875,31599,30376,29822,29691,28765,
#'  28695,28917,28203,29209,30316,29062,26977,26577,
#'  27727,28599,30513,31774,32347,34093,33736,32085,
#'  30807,28279,26873,25612,23503,22207,21388,20122,
#'  18014,15626,15006,14158,11195,7931,7640,9053,
#'  13276,17226,18918,17697,18424,17723,16706,14410,
#'  13342,14787,15183,15727,16045,14777,14267,13102,
#'  10866,9311,6933,5030,3785,3551,2848,3080,
#'  2874,2368,2681,3165,3010,3009,2721,2705,
#'  2492,2244,1971,1644,1565,1307,5027)
#'  Age <- 0:85
#'  # standard pop taken from ages 55+
#'  StPop <- c(6258,6177,6089,5995,5894,5787,5672,5552,
#'  5423,5286,5140,4985,4824,4652,4477,4293,
#'  4107,3912,3712,3502,3282,3055,2823,2591,
#'  2360,2138,1921,1710,1502,1297,1098,910,
#'  741,592,463,353,265,192,137,95,
#'  63,42,26,17,9,13) 
#'  
#'  StAge <- 55:100
#'  
#' PopExtended <- OPAG_simple(
#'  		Pop = Pop,
#'  		Age = Age,
#'  		StPop = StPop,
#'  		StAge = StAge)
#'  
#' \dontrun{
#'  plot(Age, Pop, type = 'l',xlim=c(80,100),ylim=c(0,1e4))
#' lines(0:100, PopExtended, col = "red", lty = 2)
#' }
#' stopifnot((sum(PopExtended[86:101]) - Pop[86]) == 0)

OPAG_simple <- function(Pop, Age, OAnow = max(Age), StPop, StAge, OAnew = max(StAge)){
	
	# assume single
	stopifnot(is_single(Age))
	stopifnot(is_single(StAge))
	# OAG can be less than or equal to max age
	stopifnot(OAnow %in% Age)
	# age and pop vectors must match lengths, assume ordered
	stopifnot(length(Pop) == length(Age))
	# age concordance
	#stopifnot(all(Age %in% StAge))
	
	# group pop down to OAG
	Pop      <- groupOAG(Pop, Age, OAnow)
	StPop    <- groupOAG(StPop, StAge, OAnew)
	
	# even up lengths 
	N        <- length(Pop)
	Age      <- Age[1:N]
	OAtot    <- Pop[N]
	# same for standard
	StN      <- length(StPop)
	StAge    <- StAge[1:StN]
	
	# make stadnard distribution.
	standard <- rescale.vector(StPop[StAge >= OAnow], scale = 1)
	# redistribute OAG
	PopUpper <- OAtot * standard
	# keep lower ages of Pop
	PopLower <- Pop[1:(N-1)]
	
	# graft, name, and return
	out      <- c(PopLower, PopUpper)
	Ageout   <- sort(unique(c(Age,StAge)))
	names(out) <- Ageout
	
	out
}
	











## Author: sean
## Edited 9-Dec-2017, Tim Riffe
################################################################################
##
#
##' Distributes the population of an open-ended age group at age 65 into 5-year age groups with open-ended age group at age 80.
##' @description This follows the processes set forward in the OPAG PAS spreadsheet to distribute an open-ended age group through age 65 out to an open-ended age group at age 80.
##' 
##' @param Value   numeric. A vector of demographic population counts organized into 5-year age groups with one column for s and one column for fes.
##' @param Age numeric vector of ages corresponding to the lower integer bound of the 5-year age ranges.
##' @param BirthLE numeric. A value for the life expectancy at birth for the population with one column for s and one for fes.
##' @param LTCode  string. A description of the type of life table to use for the distribution. Default "Coale-Demeny West".
##' @param LxLT numeric. A vector of Lx values from a life table that describes the population distribution at the desired age intervals with one column for s and one for fes. This must be structured in 5-year age groups up to age 80. Default is \code{NULL}.
##' @param FinalLE numeric. A vector of LE values for the final age group with one column for s and one for fes. Default is \code{NULL}.
#
##' @details 
##' @return a vector of the age distribution of the smoothed population in 5 year intervals up to and including the new open-ended age group at age 80.
##' 
##' @export
##' @author Sean Fennel
#
#opag <- function(Value, 
#		Age, 
#		e0, 
#		LTCode = "coale-demeny west", 
#		LxLT = NULL, 
#		FinalLE = NULL, 
#		Sex = "M"){
#	LTCode       <- tolower(LTCode)
#  #Set up the Lx values of the life table.
#  if(LTCode == "coale-demeny west"){
#    #Look up the referenced life table.
#    modelLT      <- getModelLT(LTCode, Sex)
#	ex           <- modelLT$ex
#    maxlevel     <- ncol(ex)
#	nLx          <- modelLT$nLx
#    #Find the largest level of the life table that is less than the given LE.
#    LevelLow <- max(which(ex[,1] < e0))
#   
#    #The OPAG spreadsheet only uses every other column (the odd columns) of Coale-Demeny
#    if(LevelLow %% 2 == 0){
#      LevelLow      <- LevelLow - 1
#    }
#    
#    LevelHigh       <- LevelLow + 2
#
#    #Interpolate between the Lx and e0 values for the levels
#    LxLTLow  <- nLx[LevelLow,]*100000
#    LxLTHigh <- nLx[LevelHigh,]*100000
#
#    E0Low           <- ex[LevelLow,1]
#
#    E0High          <- ex[LevelHigh,1]
#
#    FactorLow       <- (E0High - e0) / (E0High - E0Low)
#    FactorHigh      <- 1 - FactorLow
#
#    LxLTInit <- LxLTLow*FactorLow[1,] + LxLTHigh*FactorHigh[1,]
#    E80             <- ex[LevelLow, 18]*FactorLow + ex[LevelHigh, 18]*FactorHigh
#
#    LxLTInit2 <- c(LxLTInit[1:(length(Value[,1]))], sum(LxLTInit[(length(Value[,1])+1):length(LxLTInit)]))
#    
#    LxLTInit2 <- convertSplitTo5Year(LxLTInit)
#
#    
#    LxLT <- c(LxLTInit2[1:16,], sum(LxLTInit2[17:length(LxLTInit2[,1]),1]))
#
#    #Bring things together
#    LT5Lx <- data.frame(LxLT)
#    FinalLE      <- data.frame(E80)
#
#  }
#  else if(LTCode == "Empirical"){
#    LT5Lx      <- LxLT
#  }
#  
#  #set up different life tables for the calculations
#  LTLxPlus5    <- LT5Lx[-1,]
#  LT5LxShort   <- LT5Lx[0:(length(LT5Lx[,1])-1),]
#  LT10Lx       <- LT5LxShort + LTLxPlus5
#  LT10Lx       <- LT10Lx[c(FALSE, TRUE),]
#  
#  #set up the different age structures for the calculations below (xPlus10, x)
#  Value1Px        <- splitToSingleAges(Value[,1], Age)
#  Value10Px       <- groupAges(Value1Px, Age = seq(0,max(Age)+4, by = 1), N=10, shiftdown = 5)
#  Value10Px           <- Value10Px
#  
#  #Estimate the intrinsic growth rate.
#  r <- (log(Value10Px[6] / LT10Lx[6]) - log(Value10Px[5] / LT10Lx[5])) / 10
#  
#  #Preliminary 5-year age distribution
#  BaseValue45      <- Value[10,]
#  BaseLxValue45    <- LT5Lx[10,]
#  AgeOffset        <- seq(0, 80, by = 5) - 45
#  PrelimValues     <- (BaseValue45[, 1])*(LT5Lx[, 1]/BaseLxValue45[, 1])*exp(r[, 1] * AgeOffset)
#
#  #Update the preliminary open-ended age group.
#  PrelimValues[length(PrelimValues)] <- ((BaseValue45[, 1]) * (LT5Lx[length(LT5Lx[, 1]), 1] / BaseLxValue45[,1])* exp(r[, 1] * (AgeOffset[length(AgeOffset)] + 0.6 * E80 + 0.92)))[[1]]
#
#  
#  #Distribute preliminary based on proportion of actual totals.
#  OpenEnd <- round(PrelimValues[length(Value[,1]):(length(PrelimValues))]*(Value[length(Value[,1]),1] / sum(PrelimValues[length(Value[,1]):(length(PrelimValues))])))
#
#  #Create final data frame
#  FinalOpenEnd <- data.frame(OpenEnd)
#
#  OutputFrame <- rbind(Value[0:(length(Value[,1])-1),], FinalOpenEnd)
#  
#  return(OutputFrame)
#}