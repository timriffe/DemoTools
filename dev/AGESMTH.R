# Author: sean
# modified by TR 18 July, 2018
###############################################################################
#
#MalePop      <- c(642367, 515520, 357831, 275542, 268336, 278601, 242515, 198231, 165937, 122756, 96775, 59307, 63467, 32377, 29796, 16183, 34729)
#FemalePop    <- c(654258, 503070, 323460, 265534, 322576, 306329, 245883, 179182, 145572, 95590, 81715, 48412, 54572, 28581, 26733, 14778, 30300)
#underOnePop <- c(321476, 332782, 503070, 323460, 265534, 322576, 306329, 245883, 179182, 145572, 95590, 81715, 48412, 54572, 28581, 26733, 14778, 30300)

#' Smooth the age distribution of a population using one of five methods: Carrier-Farrag, Karup-King-Newton, Arriaga, United Nations, or Strong.
#' @description Uses the smoothing method defined by \code{Method} to generate a smoothed five year age distribution This comes from the PAS spreadsheet called AGESMTH.

#' @param Value   numeric. A vector of demographic population counts.
#' @param Age numeric vector of ages corresponding to the lower integer bound of the age range.
#' @param Method string defining the type of method. Options are "Carrier-Farrag", "KKN", "Arriaga", "UN", and "Strong".
#' @param SplitU5 boolean value indicating whether the under five age group is split into 0 & 1-4.

#' @details The Carrier-Farrag, Karup-King-Newton, and Arriaga methods do not modify the totals in each 10-year age group; the United Nations and Strong methods do. The age group structure of the output is five year age groups. The under 10 age group and the final 10-year age group are included in the output but are unable to be smoothed because of the lack of a lower or higher (respectively) 10-year interval for each of the methods.

#' @return a vector of the age distribution of the smoothed population in 5 year intervals.
#' 
#' @export
#' 
#' @examples 
#' MalePop      <- c(642367, 515520, 357831, 275542, 268336, 278601, 242515, 
#'   198231, 165937, 122756, 96775, 59307, 63467, 32377, 29796, 16183, 34729)
#' FemalePop    <- c(654258, 503070, 323460, 265534, 322576, 306329, 245883, 
#'   179182, 145572, 95590, 81715, 48412, 54572, 28581, 26733, 14778, 30300)
#' underOnePop  <- c(321476, 332782, 503070, 323460, 265534, 322576, 306329, 
#'   245883, 179182, 145572, 95590, 81715, 48412, 54572, 28581, 26733, 14778, 30300)
#' Ages <- seq(0, 80, by = 5)
#' Ages2 <- c(0,1,seq(5, 80, by = 5))
#' 
#' popAgeSmth(MalePop, Ages, "Carrier-Farrag")
#' popAgeSmth(FemalePop, Ages, "KKN")
#' popAgeSmth(underOnePop, Ages2, "Arriaga", SplitU5 = TRUE)
#' popAgeSmth(MalePop, Ages, "UN", SplitU5 = FALSE)
#' popAgeSmth(underOnePop, Ages2, Method="Strong", SplitU5 = TRUE)

# TR: SplitU5 can be superceded by detection from Age vector.
# see is.abridged()
popAgeSmth <- function(
		Value, 
		Age, 
		Method, 
		SplitU5 = FALSE){
  if (SplitU5){
    intermediateAgg <- convertSplitTo5Year(Value) #Consolidate under split under 5 group
    newAges <- seq(0, max(Age), by=5) #Create new age groups
    aggValue <- splitToSingleAges(Value, newAges) #split into single age groups
  }
  else{
    aggValue <- splitToSingleAges(Value, Age) #split into single age groups
    newAges <- seq(0, length(aggValue)-1) #Create new age groups
  }
  
  # Aggregate groups and create method inputs
  if (Method == "Carrier-Farrag" || Method == "KKN" || Method == "Arriaga" || Method == "Strong"){
    Value10PxMinus10 <- groupAges(aggValue, Age = newAges, N = 10)
    Value10Px        <- Value10PxMinus10[-1]
    Value10PxPlus10  <- Value10Px[-1]
    
    Value10PxPlus10  <- Value10PxPlus10[0:(length(Value10PxPlus10)-1)]
    Value10Px        <- Value10Px[0:(length(Value10Px)-2)]
    Value10PxMinus10 <- Value10PxMinus10[0:(length(Value10PxMinus10)-3)]
  }
  else if (Method == "UN"){
    Value5PxMinus10  <- groupAges(aggValue, Age = newAges, N = 5) #aggregate to 5 year groups
    Value5PxMinus5   <- Value5PxMinus10[-1]
    Value5Px         <- Value5PxMinus5[-1]
    Value5PxPlus5    <- Value5Px[-1]
    Value5PxPlus10   <- Value5PxPlus5[-1]
    
    Value5PxPlus10   <- Value5PxPlus10[0:(length(Value5PxPlus10)-1)]
    Value5PxPlus5    <- Value5PxPlus5[0:(length(Value5PxPlus5)-2)]
    Value5Px         <- Value5Px[0:(length(Value5Px)-3)]
    Value5PxMinus5   <- Value5PxMinus5[0:(length(Value5PxMinus5)-4)]
    Value5PxMinus10  <- Value5PxMinus10[0:(length(Value5PxMinus10)-5)]
  }
  
  # Apply method-specific formula
  if (Method == "Carrier-Farrag"){
    Value5PxPlus5    <- Value10Px / (1 + (Value10PxMinus10 / Value10PxPlus10)^(1/4))
    Value5Px         <- Value10Px - Value5PxPlus5
  }
  else if (Method == "KKN"){
    Value5Px         <- (1/2) * Value10Px + (1/16)*(Value10PxMinus10 - Value10PxPlus10)
    Value5PxPlus5    <- Value10Px - Value5Px
  }
  else if (Method == "Arriaga"){
    Value5PxPlus5    <- (-Value10PxMinus10 + 11*Value10Px + 2*Value10PxPlus10)/24
    Value5Px         <- Value10Px - Value5PxPlus5
  }
  else if (Method == "UN"){
    Value5PxPrime    <- (1/16 )* (-Value5PxMinus10 + 4*Value5PxMinus5 + 10*Value5Px + 
				                  4*Value5PxPlus5 - Value5PxMinus10)
  }
  else if (Method == "Strong"){
    Value10PxPrime   <- (Value10PxMinus10 + 2*Value10Px + Value10PxPlus10)/4
  }
  
  #Create output
  if (Method == "Carrier-Farrag" || Method == "KKN" || Method == "Arriaga"){
    #Combine sequences and create full 5-year age group structure
    SmthPop                                    <- rep(0, length(aggValue)/5)
    SmthPop[seq(3, length(SmthPop) - 4, by=2)] <- Value5Px
    SmthPop[seq(4, length(SmthPop) - 3, by=2)] <- Value5PxPlus5
    
    #Add on the under 10 and final age groups
    Value5YearGroups                           <- groupAges(aggValue, Age = newAges)
    SmthPop[seq(1, 2)]                         <- Value5YearGroups[1:2]
    SmthPop[seq(length(SmthPop)-2, length(SmthPop))] <- Value5YearGroups[(length(Value5YearGroups)-2):length(Value5YearGroups)]
  }
  else if (Method == "UN"){
    #Combine sequences and create full 5-year age group structure
    SmthPop                                    <- rep(0, length(aggValue)/5)
    SmthPop[seq(3, length(SmthPop) - 3)]       <- Value5PxPrime
    
    #Add on the under 10 and final age groups
    Value5YearGroups                           <- groupAges(aggValue, Age = newAges)
    SmthPop[seq(1, 2)]                         <- Value5YearGroups[1:2]
    SmthPop[seq(length(SmthPop)-2, length(SmthPop))] <- Value5YearGroups[(length(Value5YearGroups)-2):length(Value5YearGroups)]
  }
  else if (Method == "Strong"){
    #Combine sequences and create full 5-year age group structure
    SmthPop                                    <- rep(0, length(aggValue)/5)
    SmthPop[seq(3, length(SmthPop) - 4, by=2)] <- Value10PxPrime/2
    SmthPop[seq(4, length(SmthPop) - 3, by=2)] <- Value10PxPrime/2
    
    #Add on the under 10 and final age groups
    Value5YearGroups                           <- groupAges(aggValue, Age = newAges)
    SmthPop[seq(1, 2)]                         <- Value5YearGroups[1:2]
    SmthPop[seq(length(SmthPop)-2, length(SmthPop))] <- Value5YearGroups[(length(Value5YearGroups)-2):length(Value5YearGroups)]
  }
  
  return(SmthPop)
}