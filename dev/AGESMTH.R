# Author: sean
# modified by TR 18 July, 2018
# TODO: see individual TODO notes for functions
###############################################################################

# code simplified to enable testing. Individual methods now in /R since their examples work
Value      <- c(642367, 515520, 357831, 275542, 268336, 278601, 242515, 
		198231, 165937, 122756, 96775, 59307, 63467, 32377, 29796, 16183, 34729)
Age      <- seq(0, 80, by = 5)

smth_all <- function(Value, Age, OAG){
	data.frame(Age,Value,
			cf = carrier_farrag_smth(Value,Age,OAG),
			strong = strong_smth(Value,Age,OAG),
			arriaga = arriaga_smth(Value,Age,OAG),
			kkn = kkn_smth(Value,Age,OAG),
			un = united_nations_smth(Value,Age,OAG))
}

smth_all(Value,Age,TRUE)
smth_all(groupOAG(Value,Age,OAnew=75),seq(0,75,by=5),TRUE)

#' Smooth population in 5-year age groups using various methods
#' 
#' @description Smooth population counts in 5-year age groups using the Carrier-Farrag, 
#' Karup-King-Newton, Arriaga, United Nations, or Stong methods. Allows for imputation 
#' of values in the youngest and oldest age groups for the Carrier-Farrag, Karup-King-Newton, 
#' and United Nations methods.

#' @details The Carrier-Farrag, Karup-King-Newton, and Arriaga methods do not modify the totals 
#' in each 10-year age group, whereas the United Nations and Strong methods do. The age intervals 
#' of input data could be any integer structure (single, abridged, 5-year, other), but 
#' output is always in 5-year age groups. All methods except the United Nations methods
#' operate based on 10-year age group totals, excluding the open age group. 
#' 
#' The Carrier-Farrag, Karup-King-Newton, and United Nations methods do not produce estimates 
#' for the first and final 10-year age groups. By default, these are imputed with the original 5-year age group totals, but
#' you can also specify to impute with \code{NA}, or the results of the Arriaga or
#' Strong methods. If the terminal digit of the open age group is 5, then the terminal 10-year 
#' age group shifts down, so imputations may affect more ages in this case. Imputation can follow 
#' different methods for young and old ages.
#' 
#' Method names are simplified using \code{simplify.text} and checked against a set of plausible matches 
#' designed to give some flexibility in case you're not sure 
#' 
#' In accordance with the description of these methods in Arriaga (1994), it is advised to 
#' compare the resutls from a variety of methods. 
#' 
#' @param Value numeric vector of counts in single, abridged, or 5-year age groups.
#' @param Age integer vector of ages corresponding to the lower integer bound of the counts.
#' @param method
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}. 
#' @param ageMin integer. The lowest age included included in intermediate adjustment. Default 10. Only relevant for Strong method.
#' @param ageMax integer. The highest age class included in intermediate adjustment. Default 65. Only relevant for Strong method.
#' @param young.tail \code{NA} or character. Method to use for ages 0-9. Default \code{"original}.
#' @param old.tail \code{NA} or character. Method to use for the final age groups. Default \code{"original"}.
#' @return numeric vector of smoothed counts in 5-year age groups.
#' @export
#' 
#' @examples
#' MalePop      <- c(642367, 515520, 357831, 275542, 268336, 278601, 242515, 
#' 		198231, 165937, 122756, 96775, 59307, 63467, 32377, 29796, 16183, 34729)
#' Ages         <- seq(0, 80, by = 5)
#'
#' # names a bit flexible:
#' cf <- agesmth(Value = MalePop, 
#'		Age = Ages, 
#'		method = "Carrier-Farrag", 
#'		OAG = TRUE)
#'# old.tail be default same as young.tail
#'# "cf" also works
#'
#'# no need to specify tails for Arriaga or Strong
#'arr <- agesmth(Value = MalePop, 
#'		Age = Ages, 
#'		method = "Arriaga", 
#'		OAG = TRUE)
#'strong <- agesmth(Value = MalePop, 
#'		Age = Ages, 
#'		method = "Strong", 
#'		OAG = TRUE)
#'# other methods:
#'un <- agesmth(Value = MalePop, 
#'		Age = Ages, 
#'		method = "United Nations", 
#'		OAG = TRUE)
#'kkn <- agesmth(Value = MalePop, 
#'		Age = Ages, 
#'		method = "Karup-King-Newton", 
#'		OAG = TRUE)
#'
#'\dontrun{
#'	plot(Ages,MalePop,pch=16)
#'	lines(Ages, cf)
#'	lines(Ages, arr, col = "red")
#'	lines(Ages, strong, col = "#FF000080", lwd = 3)
#'	lines(Ages, kkn, col = "blue")
#'	lines(Ages, un, col = "magenta")
#'	legend("topright",
#'			pch=c(16,NA,NA,NA,NA,NA),
#'			lty = c(NA,1,1,1,1,1),
#'			lwd = c(NA,1,1,3,1,1),
#'			col = c("black","black","red","#FF000080","blue","magenta"),
#'			legend = c("orig 5","Carrier-Farrag","Arriaga","Strong","KKN","UN"))
#'}
#' # an extreme case:
#'  Value <- c(80626,95823,104315,115813,100796,105086,97266,116328,
#'  		75984,89982,95525,56973,78767,65672,53438,85014,
#'  		47600,64363,42195,42262,73221,30080,34391,29072,
#'  		20531,66171,24029,44227,24128,23599,82088,16454,
#'  		22628,17108,12531,57325,17220,28425,16206,17532,
#'  		65976,11593,15828,13541,8133,44696,11165,18543,
#'  		12614,12041,55798,9324,10772,10453,6773,28358,
#'  		9916,13348,8039,7583,42470,5288,5317,6582,
#'  		3361,17949,3650,5873,3279,3336,27368,1965,
#'  		2495,2319,1335,12022,1401,1668,1360,1185,
#'  		9167,424,568,462,282,6206,343,409,333,291,4137,133,169,157,89,2068,68,81,66,57)
#'  Age <- 0:99
#'  
#'  V5 <- groupAges(Value, Age=Age)
#'  Age5 <- as.integer(names(V5))
#'  cf2 <- agesmth(Value = Value, 
#'		  Age = Age, 
#'		  method = "Carrier-Farrag", 
#'		  OAG = TRUE)
#'  st2 <- agesmth(Value = Value, 
#'		  Age = Age, 
#'		  method = "Strong", 
#'		  OAG = TRUE)
#'  \dontrun{
#'  plot(Age,Value,pch=16)
#'  lines(Age,splitUniform(V5,Age=Age5,OAG=FALSE), lty=2, lwd = 2)
#'  lines(Age,splitUniform(cf2,Age=Age5,OAG=FALSE),col="blue")
#'  lines(Age,splitUniform(st2,Age=Age5,OAG=FALSE),col="red")
#'  legend("topright",
#'		  pch=c(16,NA,NA,NA),
#'		  lty=c(NA,2,1,1),
#'		  col=c("black","black","blue","red"),
#'		  lwd=c(NA,2,1,1),
#'		  legend=c("orig single","orig 5","Carrier-Farrag","Strong"))
#'}
#'  
#'# it might make sense to do this level of smoothing as intermediate step
#'# in Sprague-like situation. Compare:
#'spr1 <- sprague(Value, Age=Age,OAG=FALSE)
#'spr2 <- sprague(cf2, Age=Age5,OAG=FALSE)
#'spr3 <- sprague(st2, Age=Age5,OAG=FALSE)
#'\dontrun{
#'plot(Age,Value,pch=16, main = "Smoothing as pre-step to graduation")
#'lines(Age,spr1,lty=2)
#'lines(Age,spr2,col="blue")
#'lines(Age,spr3,col="red")
#'legend("topright",
#'		pch=c(16,NA,NA,NA),
#'		lty=c(NA,2,1,1),
#'		col=c("black","black","blue","red"),
#'		lwd=c(NA,2,1,1),
#'		legend=c("orig single","orig->Sprague","C-F->Sprague","Strong->Sprague"))
#'}

agesmth <- function(Value, 
		Age, 
		method = c("Carrier-Farrag","KKN","Arriaga","United Nations","Strong")[1], 
		OAG = TRUE, 
		ageMin = 10, ageMax = 65,
		young.tail = c("Original","Arriaga","Strong",NA)[1], 
		old.tail = young.tail){
	method     <- simplify.text(method)
	young.tail <- simplify.text(young.tail)
	old.tail   <- simplify.text(old.tail)
	
	if (missing(Age)){
		Age <- as.integer(names(Value))
	}
	stopifnot(length(Age) == length(Value))
	# carrierfarrag or cf
	if (method %in% c("cf", "carrierfarrag")){
		out <- carrier_farrag_smth(Value = Value, Age = Age, OAG = OAG)
	}
	
	# stong
	if (method == "strong"){
		out <- strong_smth(Value = Value, Age = Age, OAG = OAG, ageMin = ageMin, ageMax = ageMax)
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
	
	# -------------------------------
	# clean tails
	nas     <- is.na(out)
	if (any(nas) & (!is.na(old.tail) | !is.na(young.tail))){
		nrle         <- rle(as.integer(nas))
		original     <- groupAges(Value, Age = Age, N = 5)
		arriaga      <- arriaga_smth(Value, Age = Age, OAG = OAG)
		strong       <- strong_smth(Value, Age = Age, OAG = OAG)
		# are the final entries NAs?
		if (nrle$values[length(nrle$values)] == 1 & !is.na(old.tail)){
			nrle$values[1] <- 0
			old.ind        <- as.logical(rep(nrle$values, times = nrle$lengths))
			# do we want original values?
		    if (old.tail %in% c("o","orig","original")){
				stopifnot(length(original) == length(out))
				out[old.ind] <- original[old.ind]
			}
			# or arriaga?
		    if (old.tail == "arriaga"){
			    stopifnot(length(arriaga) == length(out))
			    out[old.ind] <- arriaga[old.ind]
		    }
			# or strong?
			if (old.tail == "strong"){
				stopifnot(length(strong) == length(out))
				out[old.ind] <- strong[old.ind]
			}
			
		}
		nrle             <- rle(as.integer(nas))
		# take care of young tail
		if (nrle$values[1] == 1 & !is.na(young.tail)){
			nrle$values[length(nrle$values)] <- 0
			young.ind        <- as.logical(rep(nrle$values, times = nrle$lengths))
			
			if (young.tail %in% c("o","orig","original")){
				stopifnot(length(original) == length(out))
				out[young.ind] <- original[young.ind]
			}
			# or arriaga?
			if (young.tail == "arriaga"){
				stopifnot(length(arriaga) == length(out))
				out[young.ind] <- arriaga[young.ind]
			}
			# or strong?
			if (young.tail == "strong"){
				stopifnot(length(strong) == length(out))
				out[young.ind] <- strong[young.ind]
			}
		}
	} # end tail operations
	
	out
}



# TR: Sean's code, slightly modified below.
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
##  if ( is.abridged(Age) | is_single(Age) ){
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