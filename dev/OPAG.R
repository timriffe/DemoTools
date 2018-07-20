# Author: sean
###############################################################################
#

#' Distributes the population of an open-ended age group at age 65 into 5-year age groups with open-ended age group at age 80.
#' @description This follows the processes set forward in the OPAG PAS spreadsheet to distribute an open-ended age group through age 65 out to an open-ended age group at age 80.
#' 
#' @param Value   numeric. A vector of demographic population counts organized into 5-year age groups with one column for males and one column for females.
#' @param Age numeric vector of ages corresponding to the lower integer bound of the 5-year age ranges.
#' @param BirthLE numeric. A value for the life expectancy at birth for the population with one column for males and one for females.
#' @param LTCode  string. A description of the type of life table to use for the distribution. Default "Coale-Demeny West".
#' @param LxLifeTable numeric. A vector of Lx values from a life table that describes the population distribution at the desired age intervals with one column for males and one for females. This must be structured in 5-year age groups up to age 80. Default is \code{NULL}.
#' @param FinalLE numeric. A vector of LE values for the final age group with one column for males and one for females. Default is \code{NULL}.

#' @details 
#' @return a vector of the age distribution of the smoothed population in 5 year intervals up to and including the new open-ended age group at age 80.
#' 
#' @export
#' 
#' @examples 

opag <- function(Value, Age, BirthLE, LTCode = "coale-demeny west", LxLifeTable = NULL, FinalLE = NULL){
  
  #Set up the Lx values of the life table.
  if(LTCode == "coale-demeny west"){
    #Look up the referenced life table.
    maleLT <- getModelLifeTable(LTCode, "M")
    femaleLT <- getModelLifeTable(LTCode, "F")
    
    #Find the largest level of the life table that is less than the given LE.
    maleLevelLow <- max(which(maleLT$ex[,1]<BirthLE[,1]))
    femaleLevelLow <- max(which(femaleLT$ex[,1]<BirthLE[,2]))
    
    #The OPAG spreadsheet only uses every other column (the odd columns) of Coale-Demeny
    if(maleLevelLow %% 2 == 0){
      maleLevelLow <- maleLevelLow - 1
    }
    
    if(femaleLevelLow %% 2 == 0){
      femaleLevelLow <- femaleLevelLow - 1
    }
    
    maleLevelHigh <- maleLevelLow + 2
    femaleLevelHigh <- femaleLevelLow + 2
    
    #Interpolate between the Lx and e0 values for the levels
    LxMaleLifeTableLow <- maleLT$nLx[maleLevelLow,]*100000
    LxFemaleLifeTableLow <- femaleLT$nLx[femaleLevelLow,]*100000
    LxMaleLifeTableHigh <- maleLT$nLx[maleLevelHigh,]*100000
    LxFemaleLifeTableHigh <- femaleLT$nLx[femaleLevelHigh,]*100000
    
    maleE0Low <- maleLT$ex[maleLevelLow,1]
    femaleE0Low <- femaleLT$ex[femaleLevelLow,1]
    maleE0High <- maleLT$ex[maleLevelHigh,1]
    femaleE0High <- femaleLT$ex[femaleLevelHigh,1]
    
    maleFactorLow <- (maleE0High - BirthLE[1]) / (maleE0High - maleE0Low)
    maleFactorHigh <- 1 - maleFactorLow
    femaleFactorLow <- (femaleE0High - BirthLE[2]) / (femaleE0High - femaleE0Low)
    femaleFactorHigh <- 1 - femaleFactorLow
    
    LxMaleLifeTableInit <- LxMaleLifeTableLow*maleFactorLow[1,] + LxMaleLifeTableHigh*maleFactorHigh[1,]
    LxFemaleLifeTableInit <- LxFemaleLifeTableLow*femaleFactorLow[1,] + LxFemaleLifeTableHigh*femaleFactorHigh[1,]
    maleE80 <- maleLT$ex[maleLevelLow, 18]*maleFactorLow + maleLT$ex[maleLevelHigh, 18]*maleFactorHigh
    femaleE80 <- femaleLT$ex[femaleLevelLow, 18]*femaleFactorLow + maleLT$ex[femaleLevelHigh, 18]*femaleFactorHigh
    
    LxMaleLifeTableInit2 <- c(LxMaleLifeTableInit[1:(length(Value[,1]))], sum(LxMaleLifeTableInit[(length(Value[,1])+1):length(LxMaleLifeTableInit)]))
    LxFemaleLifeTableInit2 <- c(LxFemaleLifeTableInit[1:(length(Value[,1]))], sum(LxFemaleLifeTableInit[(length(Value[,1])+1):length(LxFemaleLifeTableInit)]))
    
	# groupAges()
#    LxMaleLifeTableInit2 <- convertSplitTo5Year(LxMaleLifeTableInit)
#    LxFemaleLifeTableInit2 <- convertSplitTo5Year(LxFemaleLifeTableInit)
    LxMaleLifeTableInit2   <- groupAges(LxMaleLifeTableInit)
    LxFemaleLifeTableInit2 <- groupAges(LxFemaleLifeTableInit)
	
    LxMaleLifeTable <- c(LxMaleLifeTableInit2[1:16,], sum(LxMaleLifeTableInit2[17:length(LxMaleLifeTableInit2[,1]),1]))
    LxFemaleLifeTable <- c(LxFemaleLifeTableInit2[1:16,], sum(LxFemaleLifeTableInit2[17:length(LxFemaleLifeTableInit2[,1]),1]))
    
    #Bring things together
    LifeTable5Lx <- data.frame(LxMaleLifeTable, LxFemaleLifeTable)
    colnames(LifeTable5Lx) <- c("male", "female")
    
    FinalLE <- data.frame(maleE80, femaleE80)
    colnames(FinalLE) <- c("male", "female")
  }
  else if(LTCode == "Empirical"){
    LifeTable5Lx <- LxLifeTable
    colnames(LifeTable5Lx) <- c("male", "female")
    colnames(FinalLE) <- c("male", "female")
  }
  
  #set up different life tables for the calculations
  LifeTableLxPlus5 <- LifeTable5Lx[-1,]
  LifeTable5LxShort <- LifeTable5Lx[0:(length(LifeTable5Lx[,1])-1),]
  LifeTable10Lx <- LifeTable5LxShort + LifeTableLxPlus5
  LifeTable10Lx <- LifeTable10Lx[c(FALSE, TRUE),]
  
  #set up the different age structures for the calculations below (xPlus10, x)
  ValueMale1Px <- splitToSingleAges(Value[,1], Age)
  ValueFemale1Px <- splitToSingleAges(Value[,2], Age)
  ValueMale10Px <- groupAges(ValueMale1Px, Age = seq(0,max(Age)+4, by = 1), N=10, shiftdown = 5)
  ValueFemale10Px <- groupAges(ValueFemale1Px, Age = seq(0,max(Age)+4, by = 1), N=10, shiftdown = 5)
  Value10Px <- data.frame(ValueMale10Px, ValueFemale10Px)[-1,]
  colnames(Value10Px) <- c("Male", "Female")
  
  #Estimate the intrinsic growth rate.
  r <- (log(Value10Px[6,] / LifeTable10Lx[6,]) - log(Value10Px[5,] / LifeTable10Lx[5,])) / 10
  
  #Preliminary 5-year age distribution
  BaseValue45 <- Value[10,]
  BaseLxValue45 <- LifeTable5Lx[10,]
  AgeOffset <- seq(0, 80, by = 5) - 45
  PrelimValuesMale <- (BaseValue45[,1])*(LifeTable5Lx[,1]/BaseLxValue45[,1])*exp(r[,1]*AgeOffset)
  PrelimValuesFemale <- (BaseValue45[,2])*(LifeTable5Lx[,2]/BaseLxValue45[,2])*exp(r[,2]*AgeOffset)
  
  #Update the preliminary open-ended age group.
  PrelimValuesMale[length(PrelimValuesMale)] <- ((BaseValue45[,1])*(LifeTable5Lx[length(LifeTable5Lx[,1]),1]/BaseLxValue45[,1])*exp(r[,1]*(AgeOffset[length(AgeOffset)]+0.6*maleE80+0.92)))[[1]]
  PrelimValuesFemale[length(PrelimValuesFemale)] <- ((BaseValue45[,2])*(LifeTable5Lx[length(LifeTable5Lx[,2]),2]/BaseLxValue45[,2])*exp(r[,2]*(AgeOffset[length(AgeOffset)]+0.6*femaleE80+0.92)))[[1]]
  
  #Distribute preliminary based on proportion of actual totals.
  OpenEndMale <- round(PrelimValuesMale[length(Value[,1]):(length(PrelimValuesMale))]*(Value[length(Value[,1]),1] / sum(PrelimValuesMale[length(Value[,1]):(length(PrelimValuesMale))])))
  OpenEndFemale <- round(PrelimValuesFemale[length(Value[,1]):(length(PrelimValuesFemale))]*(Value[length(Value[,1]),2] / sum(PrelimValuesFemale[length(Value[,1]):(length(PrelimValuesFemale))])))
  
  #Create final data frame
  FinalOpenEnd <- data.frame(OpenEndMale, OpenEndFemale)
  colnames(FinalOpenEnd) <- c("male", "female")
  colnames(Value) <- c("male", "female")
  OutputFrame <- rbind(Value[0:(length(Value[,1])-1),], FinalOpenEnd)
  
  return(OutputFrame)
}