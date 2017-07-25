# Author: sean
###############################################################################
#
#MalePop      <- c(642367, 515520, 357831, 275542, 268336, 278601, 242515, 198231, 165937, 122756, 96775, 59307, 63467, 32377, 29796, 16183, 34729)
#FemalePop    <- c(654258, 503070, 323460, 265534, 322576, 306329, 245883, 179182, 145572, 95590, 81715, 48412, 54572, 28581, 26733, 14778, 30300)
#underOnePop <- c(321476, 332782, 503070, 323460, 265534, 322576, 306329, 245883, 179182, 145572, 95590, 81715, 48412, 54572, 28581, 26733, 14778, 30300)

#' Smooth the age distribution of a population using the Carrier-Farrag method..
#' @description Uses the Carrier-Farrag smoothing method. .... generates age-specific and summary measures; use it to evaluate the smoothing method. This comes from the PAS spreadsheet called AGESMTH.

#' @param Value   numeric. A vector of demographic population counts.
#' @param ageStruct string. An indicator of the age structure of the demographic population. Default "five year".

#' @details The Carrier-Farrag method does not modify the totals in each 10-year age group. The age group structure of the output is the same as that of the input. The possible values for the age structure are "five year", "single year", and "split under five". The under 10 age group and the final 10-year age group are included in the output but are unable to be smoothed because of the lack of a lower or higher (respectively) 10-year interval needed for the Carrier-Farrag method.

#' @return a vector of the age distribution of the smoothed population in 5 year intervals.
#' 
#' @export
#' 
#' @examples 
#' MalePop      <- c(642367, 515520, 357831, 275542, 268336, 278601, 242515, 198231, 165937, 122756, 96775, 59307, 63467, 32377, 29796, 16183, 34729)
#' FemalePop    <- c(654258, 503070, 323460, 265534, 322576, 306329, 245883, 179182, 145572, 95590, 81715, 48412, 54572, 28581, 26733, 14778, 30300)
#' underOnePop <- c(321476, 332782, 503070, 323460, 265534, 322576, 306329, 245883, 179182, 145572, 95590, 81715, 48412, 54572, 28581, 26733, 14778, 30300)
#' 
#' carrierFarragSmth(MalePop)
#' carrierFarragSmth(FemalePop, ageStruct = "five year")
#' carrierFarragSmth(underOnePop, ageStruct = "split under five")

carrierFarragSmth <- function(Value, ageStruct = "five year"){
  
  #Aggregate age groups for "single year" and "split under five"
  if (ageStruct == "single year"){
    aggValue <- convertSingleTo5Year(Value)
  }
  
  if (ageStruct == "split under five"){
    aggValue <- convertSplitTo5Year(Value)
  }
  
  if (ageStruct == "five year"){
    aggValue <- Value
  }
  
  #create vectors necessary for Carrier-Farrag formula.
  shiftedValue <- aggValue[-1]
  shiftedValue[length(aggValue)] <- 0
  initialSum <- aggValue + shiftedValue
  
  Value10PxMinus10 <- initialSum[c(TRUE,FALSE)]
  Value10Px <- Value10PxMinus10[-1]
  Value10PxPlus10 <- Value10Px[-1]
  
  Value10PxPlus10 <- Value10PxPlus10[0:(length(Value10PxPlus10)-1)]
  Value10Px <- Value10Px[0:(length(Value10Px)-2)]
  Value10PxMinus10 <- Value10PxMinus10[0:(length(Value10PxMinus10)-3)]


  #Apply Carrier-Farrag formulas
  Value5PxPlus5 = round(Value10Px / (1 + (Value10PxMinus10 / Value10PxPlus10)^(1/4)))
  Value5Px = Value10Px - Value5PxPlus5
  
  #Combine sequences and create full 5-year age group structure
  SmthPop <- rep(0, length(aggValue))
  SmthPop[seq(3, length(SmthPop) - 4, by=2)] <- Value5Px
  SmthPop[seq(4, length(SmthPop) - 3, by=2)] <- Value5PxPlus5
  SmthPop[seq(1, 2)] <- aggValue[1:2]
  SmthPop[seq(length(SmthPop)-2, length(SmthPop))] <- aggValue[(length(aggValue)-2):length(aggValue)]
  
  return(SmthPop)
}