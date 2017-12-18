# Author: sean
# Edited 9-Dec-2017, Tim Riffe
###############################################################################
#

#' Distributes the population of an open-ended age group at age 65 into 5-year age groups with open-ended age group at age 80.
#' @description This follows the processes set forward in the OPAG PAS spreadsheet to distribute an open-ended age group through age 65 out to an open-ended age group at age 80.
#' 
#' @param Value   numeric. A vector of demographic population counts organized into 5-year age groups with one column for s and one column for fes.
#' @param Age numeric vector of ages corresponding to the lower integer bound of the 5-year age ranges.
#' @param BirthLE numeric. A value for the life expectancy at birth for the population with one column for s and one for fes.
#' @param LTCode  string. A description of the type of life table to use for the distribution. Default "Coale-Demeny West".
#' @param LxLT numeric. A vector of Lx values from a life table that describes the population distribution at the desired age intervals with one column for s and one for fes. This must be structured in 5-year age groups up to age 80. Default is \code{NULL}.
#' @param FinalLE numeric. A vector of LE values for the final age group with one column for s and one for fes. Default is \code{NULL}.

#' @details 
#' @return a vector of the age distribution of the smoothed population in 5 year intervals up to and including the new open-ended age group at age 80.
#' 
#' @export
#' @author Sean Fennel

opag <- function(Value, 
		Age, 
		e0, 
		LTCode = "coale-demeny west", 
		LxLT = NULL, 
		FinalLE = NULL, 
		Sex = "M"){
	LTCode       <- tolower(LTCode)
  #Set up the Lx values of the life table.
  if(LTCode == "coale-demeny west"){
    #Look up the referenced life table.
    modelLT      <- getModelLT(LTCode, Sex)
	ex           <- modelLT$ex
    maxlevel     <- ncol(ex)
	nLx          <- modelLT$nLx
    #Find the largest level of the life table that is less than the given LE.
    LevelLow <- max(which(ex[,1] < e0))
   
    #The OPAG spreadsheet only uses every other column (the odd columns) of Coale-Demeny
    if(LevelLow %% 2 == 0){
      LevelLow      <- LevelLow - 1
    }
    
    LevelHigh       <- LevelLow + 2

    #Interpolate between the Lx and e0 values for the levels
    LxLTLow  <- nLx[LevelLow,]*100000
    LxLTHigh <- nLx[LevelHigh,]*100000

    E0Low           <- ex[LevelLow,1]

    E0High          <- ex[LevelHigh,1]

    FactorLow       <- (E0High - e0) / (E0High - E0Low)
    FactorHigh      <- 1 - FactorLow

    LxLTInit <- LxLTLow*FactorLow[1,] + LxLTHigh*FactorHigh[1,]
    E80             <- ex[LevelLow, 18]*FactorLow + ex[LevelHigh, 18]*FactorHigh

    LxLTInit2 <- c(LxLTInit[1:(length(Value[,1]))], sum(LxLTInit[(length(Value[,1])+1):length(LxLTInit)]))
    
    LxLTInit2 <- convertSplitTo5Year(LxLTInit)

    
    LxLT <- c(LxLTInit2[1:16,], sum(LxLTInit2[17:length(LxLTInit2[,1]),1]))

    #Bring things together
    LT5Lx <- data.frame(LxLT)
    FinalLE      <- data.frame(E80)

  }
  else if(LTCode == "Empirical"){
    LT5Lx      <- LxLT
  }
  
  #set up different life tables for the calculations
  LTLxPlus5    <- LT5Lx[-1,]
  LT5LxShort   <- LT5Lx[0:(length(LT5Lx[,1])-1),]
  LT10Lx       <- LT5LxShort + LTLxPlus5
  LT10Lx       <- LT10Lx[c(FALSE, TRUE),]
  
  #set up the different age structures for the calculations below (xPlus10, x)
  Value1Px        <- splitToSingleAges(Value[,1], Age)
  Value10Px       <- groupAges(Value1Px, Age = seq(0,max(Age)+4, by = 1), N=10, shiftdown = 5)
  Value10Px           <- Value10Px
  
  #Estimate the intrinsic growth rate.
  r <- (log(Value10Px[6] / LT10Lx[6]) - log(Value10Px[5] / LT10Lx[5])) / 10
  
  #Preliminary 5-year age distribution
  BaseValue45      <- Value[10,]
  BaseLxValue45    <- LT5Lx[10,]
  AgeOffset        <- seq(0, 80, by = 5) - 45
  PrelimValues     <- (BaseValue45[, 1])*(LT5Lx[, 1]/BaseLxValue45[, 1])*exp(r[, 1] * AgeOffset)

  #Update the preliminary open-ended age group.
  PrelimValues[length(PrelimValues)] <- ((BaseValue45[, 1]) * (LT5Lx[length(LT5Lx[, 1]), 1] / BaseLxValue45[,1])* exp(r[, 1] * (AgeOffset[length(AgeOffset)] + 0.6 * E80 + 0.92)))[[1]]

  
  #Distribute preliminary based on proportion of actual totals.
  OpenEnd <- round(PrelimValues[length(Value[,1]):(length(PrelimValues))]*(Value[length(Value[,1]),1] / sum(PrelimValues[length(Value[,1]):(length(PrelimValues))])))

  #Create final data frame
  FinalOpenEnd <- data.frame(OpenEnd)

  OutputFrame <- rbind(Value[0:(length(Value[,1])-1),], FinalOpenEnd)
  
  return(OutputFrame)
}