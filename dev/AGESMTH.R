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
 MalePop      <- c(642367, 515520, 357831, 275542, 268336, 278601, 242515, 
   198231, 165937, 122756, 96775, 59307, 63467, 32377, 29796, 16183, 34729)

 FemalePop    <- c(654258, 503070, 323460, 265534, 322576, 306329, 245883, 
   179182, 145572, 95590, 81715, 48412, 54572, 28581, 26733, 14778, 30300)
 underOnePop  <- c(321476, 332782, 503070, 323460, 265534, 322576, 306329, 
   245883, 179182, 145572, 95590, 81715, 48412, 54572, 28581, 26733, 14778, 30300)
 Ages <- seq(0, 80, by = 5)
 Ages2 <- c(0,1,seq(5, 80, by = 5))
#' 
#' popAgeSmth(MalePop, Ages, "Carrier-Farrag")
#' popAgeSmth(FemalePop, Ages, "KKN")
#' popAgeSmth(underOnePop, Ages2, "Arriaga", SplitU5 = TRUE)
#' popAgeSmth(MalePop, Ages, "UN", SplitU5 = FALSE)
#' popAgeSmth(underOnePop, Ages2, Method="Strong", SplitU5 = TRUE)

Value <- MalePop
Age   <- Ages
OAG   <- TRUE
CFmales <- carrier_farrag_smth(MalePop, Ages, TRUE) 

# spreadsheet results
CFtest <- c(NA,NA,346290,287083,285855,261082,237937,
202809,162973,125720,88730,67352,55187,40657,NA,NA,NA)
all(round(CFmales) - CFtest == 0, na.rm = TRUE)

carrier_farrag_smth <- function(Value, 
		Age, 
		OAG = TRUE){
	
	N <- length(Value)
	if (OAG){
		Value[N] <- NA
	}
	
	Value5     <- groupAges(Value, Age = Age, N = 5)
	N          <- length(Value5)
	Age5       <- as.integer(names(Value5))
	
	# get staggered vectors
	Value5LL   <- shift.vector(Value5, -2, fill = NA)
	Value5L    <- shift.vector(Value5, -1, fill = NA)
	Value5R    <- shift.vector(Value5, 1, fill = NA)
	Value5RR   <- shift.vector(Value5, 2, fill = NA)
	Value5RRR  <- shift.vector(Value5, 3, fill = NA)
	
	# this is the funny C-F operation
	ValuePert  <-  (Value5 + Value5R) / 
			(1 + ((Value5RRR + Value5RR) / (Value5L + Value5LL))^.25)
	
	# indices for switching behavior on and off.
	# need same nr of evens and odds.
	NN         <- N
	if (N %% 2 == 1){
		NN        <- N + 1
		ValuePert <- c(ValuePert, NA)
		Value5    <- c(Value5, NA)
		Value5L   <- c(Value5L, NA)
	}
	
	inds       <- 1:NN
	odds       <- inds %% 2 == 1
	evens      <- !odds
	
	# produce results vector
	out        <- Value5 * NA 
	out[evens] <- ValuePert[evens]
	
	# make sure sum(odds) == sum(evens)
	out[odds]  <- (Value5 + Value5L)[odds] - ValuePert[evens]
	
	# cut back down (depending) and name
	out        <- out[1:N]
	names(out) <- Age5
    out
	# tail behavior will be controlled in top level function.
}
KKNtest <- c(NA,NA,354871,278502,285508,261429,236513 ,
204233,162138,126555,90094,65988,54803,41041,NA,NA,NA)

KKNmales <- kkn_smth(MalePop, Ages, TRUE)
all(round(KKNmales) - KKNtest == 0, na.rm = TRUE)

kkn_smth <- function(Value, 
		Age, 
		OAG = TRUE){
	
	N <- length(Value)
	if (OAG){
		Value[N] <- NA
	}
	
	Value5     <- groupAges(Value, Age = Age, N = 5)
	N          <- length(Value5)
	Age5       <- as.integer(names(Value5))
	
	# get staggered vectors
	Value5LLL  <- shift.vector(Value5, -3, fill = NA)
	Value5LL   <- shift.vector(Value5, -2, fill = NA)
	Value5L    <- shift.vector(Value5, -1, fill = NA)
	Value5R    <- shift.vector(Value5, 1, fill = NA)
	Value5RR   <- shift.vector(Value5, 2, fill = NA)
	Value5RRR  <- shift.vector(Value5, 3, fill = NA)
	
	# this is the funny KNN operation
	ValuePert  <-  (Value5 + Value5L) / 2 + (Value5R + Value5RR - Value5LL - Value5LLL) / 16

	# indices for switching behavior on and off.
	# need same nr of evens and odds.
	NN         <- N
	if (N %% 2 == 1){
		NN        <- N + 1
		ValuePert <- c(ValuePert, NA)
		Value5    <- c(Value5, NA)
		Value5L   <- c(Value5L, NA)
	}
	
	inds       <- 1:NN
	odds       <- inds %% 2 == 1
	evens      <- !odds
	
	# produce results vector
	out        <- Value5 * NA 
	out[odds]  <- ValuePert[odds]
	
	# make sure sum(odds) == sum(evens)
	out[evens]  <- (Value5 + Value5L)[odds] - ValuePert[odds]
	
	# cut back down (depending) and name
	out        <- out[1:N]
	names(out) <- Age5
	out
}


AMales <- arriaga_smth(MalePop, Ages, TRUE)
Atest <- c(662761, 495126, 345744, 287629, 285919, 261018, 237469, 203277, 
161733, 126960, 88586, 67496, 54587, 41257, 28790, 17189,34729 ) 
all(round(AMales) - Atest == 0, na.rm = TRUE)


# TODO: test this to vectors of different lengths: didn't work for age 70.
# need more general indexing solution.
Atest2        <- c(662761, 495126, 345744, 287629, 285919, 261018, 237469, 203277, 
		161733, 126960, 88586, 67496, 54587, 41257, 28790, 17189+34729 ) 
MalePop2      <- c(642367, 515520, 357831, 275542, 268336, 278601, 242515, 
		198231, 165937, 122756, 96775, 59307, 63467, 32377, 29796, 16183 + 34729)
AMales2 <- arriaga_smth(MalePop2, seq(0, 75, by = 5), TRUE)



arriaga_smth <- function(Value, 
		Age, 
		OAG = TRUE){
	
	N <- length(Value)
	if (OAG){
		OAGvalue <- Value[N]
		Value[N] <- NA
	}
	
	Value5     <- groupAges(Value, Age = Age, N = 5)
	N          <- length(Value5)
	Age5       <- as.integer(names(Value5))
	
	# get staggered vectors
	Value5LLLL <- shift.vector(Value5, -4, fill = NA)
	Value5LLL  <- shift.vector(Value5, -3, fill = NA)
	Value5LL   <- shift.vector(Value5, -2, fill = NA)
	Value5L    <- shift.vector(Value5, -1, fill = NA)
	Value5R    <- shift.vector(Value5, 1, fill = NA)
	Value5RR   <- shift.vector(Value5, 2, fill = NA)
	Value5RRR  <- shift.vector(Value5, 3, fill = NA)
	Value5RRRR <- shift.vector(Value5, 4, fill = NA)
	
	# thValue5RRR  <- shift.vector(Value5, 3, fill = NA)is is the funny KNN operation

	ValuePertYoung <- (8*(Value5 + Value5R) + 5 * (Value5L + Value5LL) - Value5LLL - Value5LLLL) / 24
	ValuePert      <- (-Value5RRR - Value5RR + 11 * (Value5R + Value5) + 2 * (Value5L + Value5LL)) / 24
	ValuePertOld   <- (8*(Value5 + Value5L) + 5 * (Value5R + Value5RR) - Value5RRR - Value5RRRR) / 24
	oldi           <- (N-2):N
	youngi         <- 1:3
	# cbind(ValuePertYoung, ValuePert, ValuePertOld)
	# combine into one vector.
	ValuePert[youngi] <- ValuePertYoung[youngi]
	ValuePert[oldi] <- ValuePertOld[oldi]
	
	
	# indices for switching behavior on and off.
	# need same nr of evens and odds.
	NN         <- N
	if (N %% 2 == 1){
		NN        <- N + 1
		ValuePert <- c(ValuePert, NA)
		Value5    <- c(Value5, NA)
		Value5L   <- c(Value5L, NA)
		Value5R   <- c(Value5R, NA)
	}
	
	inds       <- 1:NN
	evens      <- inds %% 2 == 0 & inds < (N-1)
	evens[N-2] <- TRUE
	
	# produce results vector
	out        <- Value5 * NA 
	out[evens] <- ValuePert[evens]
	# tougher
    odds       <- inds %% 2 == 1 & inds < N & ! evens
	out[odds]  <- (Value5 + Value5L)[odds] - ValuePert[which(odds)+1]
	# now final value might have to look 'up'
    final      <- !odds & !evens
	out[final] <- (Value5 + Value5R)[final] - shift.vector(ValuePert, 1)[final]
	# replace OAG if relevant:
    if (OAG){
		out[N]     <- OAGvalue
	}
	# cut back down (depending) and name
	out        <- out[1:N]
	names(out) <- Age5
	out
}

untest <- c(NA,NA,364491,279123,268724,272228,243638,200923,162752,126304,
91662,67432,54677,38833,NA,NA,NA)
all(round(united_nations_smth(MalePop,Ages,TRUE)) - untest == 0, na.rm = TRUE)

united_nations_smth <- function(Value, 
		Age, 
		OAG = TRUE){
	N <- length(Value)
	if (OAG){
		Value[N] <- NA
	}
	
	Value5     <- groupAges(Value, Age = Age, N = 5)
	N          <- length(Value5)
	Age5       <- as.integer(names(Value5))
	
	# get staggered vectors
	Value5LL   <- shift.vector(Value5, -2, fill = NA)
	Value5L    <- shift.vector(Value5, -1, fill = NA)
	Value5R    <- shift.vector(Value5, 1, fill = NA)
	Value5RR   <- shift.vector(Value5, 2, fill = NA)
	
	# this is the funny KNN operation
	# B11 is central
	out  <-  (-Value5RR + 4 * (Value5L + Value5R) + 10 * Value5 - Value5LL) / 16
	
	# cut back down (depending) and name
	names(out) <- Age5
	out
}

strongtest <- c(646617, 511270, 386889, 317345, 273736, 240058, 218645, 188297, 
153931, 124347, 93254, 71858, 53594, 39721, 27887, 18092,34729 
) 
# differences due to intermediate rounding in spreadsheet (bad practice IMO)
all(abs(strong_smth(MalePop,Ages,TRUE) - strongtest) < 1, na.rm = TRUE)

strong_smth <- function(Value, 
		Age, 
		OAG = TRUE,
		minA = 10,
		maxA = 65){
	
	stopifnot(maxA %% 10 == 0 & minA %% 10 == 0)
	N     <- length(Value)
	Tot   <- sum(Value)
	
	# possibly drop Open Age to something divisible by 10 just in case.	
	
	Value5     <- groupAges(Value, Age = Age, N = 5)
	
	# what OAG is a strange digit? Then take OAG after grouping.
	if (OAG){
		OAGvalue <- Value5[length(Value5)]
		Value[N] <- NA
		
	}
	Value10    <- groupAges(Value, Age = Age, N = 10)
	
	
	N          <- length(Value5)
	Age5       <- as.integer(names(Value5))
	Age10      <- as.integer(names(Value10))
	
	# subtotal
	indsub     <- Age10 >= minA & Age10 <= maxA
	SubTot     <- sum(Value10[indsub])
	
#	# get staggered vectors
	Value10L    <- shift.vector(Value10, -1, fill = NA)
	Value10R    <- shift.vector(Value10, 1, fill = NA)
	# this is the funny KNN operation
	# B11 is central
	Value10Pert <- (Value10 * 2 + Value10L + Value10R) / 4
	Value10Pert[is.na(Value10Pert)] <- Value10[is.na(Value10Pert)]
	
	# rescale ages between min and max to sum to original
	Value10Adj  <- Value10Pert
	Value10Adj[indsub] <- Value10Adj[indsub] * SubTot / sum(Value10Adj[indsub])
	
	# again get staggered vectors
	V10adjLL    <- shift.vector(Value10Adj, -2, fill = NA)
	V10adjL     <- shift.vector(Value10Adj, -1, fill = NA)
	V10adjR     <- shift.vector(Value10Adj, 1, fill = NA)
	V10adjRR    <- shift.vector(Value10Adj, 2, fill = NA)
	
	# alternating calc, with differences at tails
	V5evens     <- (-V10adjR + 11 * Value10Adj + 2 * V10adjL) / 24
	# tails different
	V5evens[1]  <- (8 * Value10Adj[1] + 5 * V10adjL[1] - V10adjLL[1]) / 24
	lastind     <- which(is.na(V5evens))[1]
	V5evens[lastind] <-  Value10Adj[lastind] - (8 * Value10Adj[lastind] + 5 * V10adjR[lastind] - V10adjRR[lastind]) / 24
	# odds are complement
	V5odds      <- Value10Adj - V5evens
	
	out         <- Value5 * NA
	names(out)  <- Age5
 	evens       <- 1:(length(V5odds)-1) * 2
	odds        <- evens - 1
	out[odds]   <- V5odds[1:lastind]
	out[evens]  <- V5evens[1:lastind]
	
	
    # keep but rescale?
    if (OAG){
		out[N] <- OAGvalue
	}
	
	# what if OAis e.g. 85?
	out
}

string <- "Carrier_Farrag"
simplifytext <- function(string){
	lower <- tolower(string)
	sub("[^[:alpha:]]+", "", lower)
}


agesmth <- function(Value, Age, method = "Carrier-Farrag", OAG = TRUE, minA = 10, maxA = 70){
	method <- simplifytext(method)
	# carrierfarrag or cf
	if (method %in% c("cf", "carrierfarrag")){
		out <- carrier_farrag_smth(Value = Value, Age = Age, OAG = OAG)
	}
	
	# stong
	if (method == "strong"){
		out <- strong_smth(Value = Value, Age = Age, OAG = OAG, minA = minA, maxA = maxA)
	}
	
	# un or unitednations
	if (method %in% c("un", "unitednations")){
		out <- united_nations_smth(Value = Value, Age = Age, OAG = OAG)
	}
	
	# arriaga
	if (method  == "arriaga"){
		out <- arriaga_smth(Value = Value, Age = Age, OAG = OAG)
	}
	
	# kkn kkingnewton karupkingnewton
	if (method %in% c("kkn", "kkingnewton", "karupkingnewton")){
		out <- kkn_smth(Value = Value, Age = Age, OAG = OAG)
	}
	
	out
}

agesmth(MalePop, Ages, "Carrier-Farrag", TRUE)
# TR: SplitU5 can be superceded by detection from Age vector.
# see is.abridged()
#popAgeSmth <- function(
#		Value, 
#		Age, 
#		Method,
#		OAG = TRUE){
#  # TR, move to single ages FROM 5-year age groups. 
#  Value5     <- groupAges(Value, Age = Age, N = 5)
#  Value10    <- groupAges(Value, Age = Age, N = 10)
#  Age5       <- as.integer(names(Value5))
#  Value1     <- splitUniform(Value5, Age = Age5, OAG = OAG)
#  
#  # vectors for "Carrier-Farrag", "KKN" , "Arriaga", "Strong"
#  Value10R   <- shift.vector(Value10, -1, fill = NA) #Value10PxMinus10
#  Value10L   <- shift.vector(Value10, 1, fill = NA)  #Value10PxPlus10
#  
#  # not sure if needed
#  Value5LL   <- shift.vector(Value5, -2, fill = NA)
#  Value5L    <- shift.vector(Value5, -1, fill = NA)
#  Value5R    <- shift.vector(Value5, 1, fill = NA)
#  Value5RR   <- shift.vector(Value5, 2, fill = NA)
#  #rbind(Value10L,Value10,Value10R)
#	
##  if ( is.abridged(Age) | is.single(Age) ){
##    #intermediateAgg <- convertSplitTo5Year(Value) #Consolidate under split under 5 group
##	intermediateAgg <- groupAges(Value,Age=Age)
##    newAges         <- as.integer(names(intermediateAgg)) 
##	# TR: by default final age group held as-is.
##    aggValue        <- splitUniform(Value, Age=newAges,OAG=OAG) #split into single age groups
##  }
##  else{
##    aggValue <- splitUniform(Value, Age = Age) #split into single age groups
##    newAges <- seq(0, length(aggValue)-1) #Create new age groups
##  }
#  
#  # Aggregate groups and create method inputs
#  if (Method %in% c("Carrier-Farrag", "KKN" , "Arriaga", "Strong")){
#    Value10PxMinus10 <- groupAges(aggValue, Age = newAges, N = 10)
#    Value10Px        <- Value10PxMinus10[-1]
#    Value10PxPlus10  <- Value10Px[-1]
#    
#    Value10PxPlus10  <- Value10PxPlus10[0:(length(Value10PxPlus10)-1)]
#    Value10Px        <- Value10Px[0:(length(Value10Px)-2)]
#    Value10PxMinus10 <- Value10PxMinus10[0:(length(Value10PxMinus10)-3)]
#  }
#  
#  # TR: this'll take a while to sort through
#  else if (Method == "UN"){
#    Value5PxMinus10  <- groupAges(aggValue, Age = newAges, N = 5) #aggregate to 5 year groups
#	# Value5
#
#	
#    Value5PxMinus5   <- Value5PxMinus10[-1]
#    Value5Px         <- Value5PxMinus5[-1]
#    Value5PxPlus5    <- Value5Px[-1]
#    Value5PxPlus10   <- Value5PxPlus5[-1]
#    
#    Value5PxPlus10   <- Value5PxPlus10[0:(length(Value5PxPlus10)-1)]
#    Value5PxPlus5    <- Value5PxPlus5[0:(length(Value5PxPlus5)-2)]
#    Value5Px         <- Value5Px[0:(length(Value5Px)-3)]
#    Value5PxMinus5   <- Value5PxMinus5[0:(length(Value5PxMinus5)-4)]
#    Value5PxMinus10  <- Value5PxMinus10[0:(length(Value5PxMinus10)-5)]
#  }
#  
#  # Apply method-specific formula
#  if (Method == "Carrier-Farrag"){
#	  
#    Value5PxPlus5    <- Value10Px / (1 + (Value10PxMinus10 / Value10PxPlus10)^(1/4))
#    Value5Px         <- Value10Px - Value5PxPlus5
#  }
#  
#  else if (Method == "KKN"){
#    Value5Px         <- (1/2) * Value10Px + (1/16)*(Value10PxMinus10 - Value10PxPlus10)
#    Value5PxPlus5    <- Value10Px - Value5Px
#  }
#  else if (Method == "Arriaga"){
#    Value5PxPlus5    <- (-Value10PxMinus10 + 11*Value10Px + 2*Value10PxPlus10)/24
#    Value5Px         <- Value10Px - Value5PxPlus5
#  }
#  else if (Method == "UN"){
#    Value5PxPrime    <- (1/16 )* (-Value5PxMinus10 + 4*Value5PxMinus5 + 10*Value5Px + 
#				                  4*Value5PxPlus5 - Value5PxMinus10)
#  }
#  else if (Method == "Strong"){
#    Value10PxPrime   <- (Value10PxMinus10 + 2*Value10Px + Value10PxPlus10)/4
#  }
#  
#  #Create output
#  if (Method %in% c("Carrier-Farrag", "KKN", "Arriaga")){
#    #Combine sequences and create full 5-year age group structure
#    SmthPop                                    <- rep(0, length(aggValue)/5)
#    SmthPop[seq(3, length(SmthPop) - 4, by=2)] <- Value5Px
#    SmthPop[seq(4, length(SmthPop) - 3, by=2)] <- Value5PxPlus5
#    
#    #Add on the under 10 and final age groups
#    Value5YearGroups                           <- groupAges(aggValue, Age = newAges)
#    SmthPop[seq(1, 2)]                         <- Value5YearGroups[1:2]
#    SmthPop[seq(length(SmthPop)-2, length(SmthPop))] <- Value5YearGroups[(length(Value5YearGroups)-2):length(Value5YearGroups)]
#  }
#  else if (Method == "UN"){
#    #Combine sequences and create full 5-year age group structure
#    SmthPop                                    <- rep(0, length(aggValue)/5)
#    SmthPop[seq(3, length(SmthPop) - 3)]       <- Value5PxPrime
#    
#    #Add on the under 10 and final age groups
#    Value5YearGroups                           <- groupAges(aggValue, Age = newAges)
#    SmthPop[seq(1, 2)]                         <- Value5YearGroups[1:2]
#    SmthPop[seq(length(SmthPop)-2, length(SmthPop))] <- Value5YearGroups[(length(Value5YearGroups)-2):length(Value5YearGroups)]
#  }
#  else if (Method == "Strong"){
#    #Combine sequences and create full 5-year age group structure
#    SmthPop                                    <- rep(0, length(aggValue)/5)
#    SmthPop[seq(3, length(SmthPop) - 4, by=2)] <- Value10PxPrime/2
#    SmthPop[seq(4, length(SmthPop) - 3, by=2)] <- Value10PxPrime/2
#    
#    #Add on the under 10 and final age groups
#    Value5YearGroups                           <- groupAges(aggValue, Age = newAges)
#	
#	# TR why in 2 steps?
#    SmthPop[seq(1, 2)]                         <- Value5YearGroups[1:2]
#    SmthPop[seq(length(SmthPop)-2, length(SmthPop))] <- Value5YearGroups[(length(Value5YearGroups)-2):length(Value5YearGroups)]
#  }
#  
#  return(SmthPop)
#}