
#' The Carrier-Farrag method of population count smoothing
#' @description Smooth population counts in 5-year age groups.
#' @details This method does not account for ages < 10 nor for the 10 year age interval prior to the open age group. These are returned imputed with \code{NA}. Age classes must be cleanly groupable to 5-year age groups. Smoothed counts are constrained to sum to original totals in 10-year age groups.
#' @param Value numeric vector of counts in single, abridged, or 5-year age groups.
#' @param Age numeric vector of ages corresponding to the lower integer bound of the counts.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}.
#' @return numeric vector of smoothed counts in 5-year age groups.
#' @export
#' @examples
#' # from PASEX AGESMTH
#' Ages         <- seq(0, 80, by = 5)
#' KKNtest <- c(NA,NA,354871,278502,285508,261429,236513 ,
#' 		204233,162138,126555,90094,65988,54803,41041,NA,NA,NA)
#'
#' CFmales <- smooth_age_5_cf(pop5m_pasex, Ages, TRUE)
#' CFtest <- c(NA,NA,346290,287083,285855,261082,237937,
#' 202809,162973,125720,88730,67352,55187,40657,NA,NA,NA)
#' all(round(CFmales) - CFtest == 0, na.rm = TRUE)
#' \dontrun{
#' plot(Ages, pop5m_pasex)
#' lines(as.integer(names(CFmales)),CFmales)
#' }
#' @references
#' \insertRef{carrier1959reduction}{DemoTools}
#' \insertRef{PAS}{DemoTools}

smooth_age_5_cf <- function(Value,
                                Age,
                                OAG = TRUE) {
  if (as.character(match.call()[[1]]) == "carrier_farrag_smth") {
    warning("please use smooth_age_5_cf() instead of carrier_farrag_smth().", call. = FALSE)
  }
  
  # these values are not used, it's just for lengths, and to make sure we
  # end on an even 10. Technically we could even provide data in 10-year
  # age groups and it'd still not break.
  Value1     <- graduate_uniform(Value = Value, Age = Age, OAG = OAG)
  Value5     <-
    groupAges(Value1, Age = as.integer(names(Value1)), N = 5)
  N          <- length(Value5)
  Age5       <- as.integer(names(Value5))
  
  # would need to move this up to ensure?
  # or in case of 85+ would we want to keep 80-84, 85+ as-is?
  Value10    <- groupAges(Value, Age = Age, N = 10)
  
  # what OAG is a strange digit? Then take OAG after grouping.
  if (OAG) {
    OAGvalue <- Value10[length(Value10)]
    Value10[length(Value10)] <- NA
  }
  
  v10R       <- shift.vector(Value10, 1, fill = NA)
  v10L       <- shift.vector(Value10,-1, fill = NA)
  vevens     <- Value10 / (1 + (v10R / v10L) ^ .25)
  vodds      <- Value10 - vevens
  
  out        <- Value5 * NA
  # cut back down (depending) and name
  interleaf  <- c(rbind(vodds, vevens))
  n          <- min(c(length(interleaf), N))
  out[1:n]   <- interleaf[1:n]
  
  out
  # tail behavior will be controlled in top level function.
}

#' @export
#' @rdname smooth_age_5_cf
carrier_farrag_smth <- smooth_age_5_cf

#' Karup-King-Newton method of population count smoothing
#' @description Smooth population counts in 5-year age groups.
#' @details This method does not account for ages < 10 nor for the 10 year age interval prior to the open age group. These are returned imputed with \code{NA}. Age classes must be cleanly groupable to 5-year age groups. Smoothed counts are constrained to sum to original totals in 10-year age groups.
#' @param Value numeric vector of counts in single, abridged, or 5-year age groups.
#' @param Age numeric vector of ages corresponding to the lower integer bound of the counts.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}.
#' @return numeric vector of smoothed counts in 5-year age groups.
#' @export
#' @examples
#' Ages         <- seq(0, 80, by = 5)
#' KKNtest <- c(NA,NA,354871,278502,285508,261429,236513 ,
#' 		204233,162138,126555,90094,65988,54803,41041,NA,NA,NA)
#'
#' KKNmales <- smooth_age_5_kkn(pop5m_pasex, Ages, TRUE)
#' all(round(KKNmales) - KKNtest == 0, na.rm = TRUE)
#' \dontrun{
#' plot(Ages, pop5m_pasex)
#' lines(as.integer(names(KKNmales)),KKNmales)
#' }
#' @references
#' \insertRef{carrier1959reduction}{DemoTools}
# plz add PASEX citation, and add above to bibtex

smooth_age_5_kkn <- function(Value,
                     Age,
                     OAG = TRUE) {
  if (as.character(match.call()[[1]]) == "kkn_smth") {
    warning("please use smooth_age_5_kkn() instead of kkn_smth().", call. = FALSE)
  }
  
  # these values are not used, it's just for lengths, and to make sure we
  # end on an even 10. Technically we could even provide data in 10-year
  # age groups and it'd still not break.
  Value1     <- graduate_uniform(Value = Value, Age = Age, OAG = OAG)
  Value5     <-
    groupAges(Value1, Age = as.integer(names(Value1)), N = 5)
  N          <- length(Value5)
  Age5       <- as.integer(names(Value5))
  
  # would need to move this up to ensure?
  # or in case of 85+ would we want to keep 80-84, 85+ as-is?
  Value10    <- groupAges(Value, Age = Age, N = 10)
  
  # what OAG is a strange digit? Then take OAG after grouping.
  if (OAG) {
    OAGvalue <- Value10[length(Value10)]
    Value10[length(Value10)] <- NA
  }
  
  v10R       <- shift.vector(Value10, 1, fill = NA)
  v10L       <- shift.vector(Value10,-1, fill = NA)
  
  # this is the KNN operation
  vodds      <- Value10 / 2 + (v10R - v10L) / 16
  # constrained in 10-year age groups
  vevens     <- Value10 - vodds
  # stagger odd even 5s
  interleaf  <- c(rbind(vodds, vevens))
  # produce results vector
  out        <- Value5 * NA
  n          <- min(c(length(interleaf), N))
  out[1:n]   <- interleaf[1:n]
  
  out
}
#' @export
#' @rdname smooth_age_5_kkn
kkn_smth <- smooth_age_5_kkn


#' E. Arriaga's method of population count smoothing
#' @description Smooth population counts in 5-year age groups.
#' @details The open age group is aggregated down to be evenly divisible by 10.
#' This method accounts for the youngest and oldest age groups. Age classes must be cleanly groupable to 5-year age groups. 10-year age groups are constrained to sum to their original totals.
#' @param Value numeric vector of counts in single, abridged, or 5-year age groups.
#' @param Age numeric vector of ages corresponding to the lower integer bound of the counts.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}.
#' @return numeric vector of smoothed counts in 5-year age groups.
#' @export
#' @examples
#' Ages         <- seq(0, 80, by = 5)
#' AMales       <- smooth_age_5_arriaga(Value = pop5m_pasex, Age = Ages, OAG = TRUE)
#' # PAS spreadsheet result:
#' Atest        <- c(662761, 495126, 345744, 287629, 285919, 261018, 237469, 203277,
#' 161733, 126960, 88586, 67496, 54587, 41257, 28790, 17189, 34729)
#' all(round(AMales) - Atest == 0, na.rm = TRUE)
#' \dontrun{
#' plot(Ages, pop5m_pasex)
#' lines(as.integer(names(AMales)),AMales)
#' }
#' @references
#' \insertRef{arriaga1994population}{DemoTools}
#' \insertRef{arriaga1968new}{DemoTools}


smooth_age_5_arriaga <- function(Value,
                         Age,
                         OAG = TRUE) {
  if (as.character(match.call()[[1]]) == "arriaga_smth") {
    warning("please use smooth_age_5_arriaga() instead of arriaga_smth().", call. = FALSE)
  }
  
  # these values are not used, it's just for lengths, and to make sure we
  # end on an even 10. Technically we could even provide data in 10-year
  # age groups and it'd still not break.
  Value1     <- graduate_uniform(Value = Value, Age = Age, OAG = OAG)
  Value5     <-
    groupAges(Value1, Age = as.integer(names(Value1)), N = 5)
  N          <- length(Value5)
  Age5       <- as.integer(names(Value5))
  
  # would need to move this up to ensure?
  # or in case of 85+ would we want to keep 80-84, 85+ as-is?
  Value10    <- groupAges(Value, Age = Age, N = 10)
  
  # what OAG is a strange digit? Then take OAG after grouping.
  if (OAG) {
    OAGvalue <- Value5[length(Value5)]
    Value10[length(Value10)] <- NA
    Value5[length(Value5)]   <- NA
  }
  
  # again get staggered vectors
  Value10LL   <- shift.vector(Value10,-2, fill = NA)
  Value10L    <- shift.vector(Value10,-1, fill = NA)
  Value10R    <- shift.vector(Value10, 1, fill = NA)
  Value10RR   <- shift.vector(Value10, 2, fill = NA)
  
  # alternating calc, with differences at tails
  vevens      <- (-Value10R + 11 * Value10 + 2 * Value10L) / 24
  # tails different
  vevens[1]   <-
    (8 * Value10[1] + 5 * Value10L[1] - Value10LL[1]) / 24
  lastind     <- which(is.na(vevens))[1]
  vevens[lastind] <-
    Value10[lastind] - (8 * Value10[lastind] + 5 * Value10R[lastind] - Value10RR[lastind]) / 24
  # odds are complement
  vodds       <- Value10 - vevens
  
  # prepare output
  interleaf  <- c(rbind(vodds, vevens))
  # produce results vector
  out        <- Value5 * NA
  n          <- min(c(length(interleaf), N))
  out[1:n]   <- interleaf[1:n]
  
  # if OA ends in 5, then we can save penultimate value too
  na.i       <- is.na(out)
  out[na.i]  <- Value5[na.i]
  if (OAG) {
    out[N] <- OAGvalue
  }
  
  out
}
#' @export
#' @rdname smooth_age_5_arriaga
arriaga_smth <- smooth_age_5_arriaga

#' The old United Nations method of population count smoothing
#' @description Smooth population counts in 5-year age groups.
#' @details The open age group is aggregated down to be evenly divisible by 10. This method accounts for the youngest and oldest age groups. Age classes must be cleanly groupable to 5-year age groups. Counts are not constrained in 10-year age groups, except 10-year young and old tails, which are unaffected.
#' @param Value numeric vector of counts in single, abridged, or 5-year age groups.
#' @param Age numeric vector of ages corresponding to the lower integer bound of the counts.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}.
#' @return numeric vector of smoothed counts in 5-year age groups.
#' @export
#' @examples
#' Ages         <- seq(0, 80, by = 5)
#' un_test <- c(NA,NA,364491,279123,268724,272228,243638,200923,162752,126304,
#' 		91662,67432,54677,38833,NA,NA,NA)
#' un_result <- smooth_age_5_un(pop5m_pasex, Ages, TRUE)
#' all(round(un_result) - un_test == 0, na.rm = TRUE)
#' \dontrun{
#' plot(Ages, pop5m_pasex)
#' lines(as.integer(names(un_result)),un_result)
#' }
#' @references
#' \insertRef{carrier1959reduction}{DemoTools}
#' \insertRef{arriaga1994population}{DemoTools}


smooth_age_5_un <- function(Value,
                                Age,
                                OAG = TRUE) {
  if (as.character(match.call()[[1]]) == "united_nations_smth") {
    warning("please use smooth_age_5_un() instead of united_nations_smth().", call. = FALSE)
  }
  
  N <- length(Value)
  if (OAG) {
    Value[N] <- NA
  }
  
  Value5     <- groupAges(Value, Age = Age, N = 5)
  N          <- length(Value5)
  Age5       <- as.integer(names(Value5))
  
  # get staggered vectors
  Value5LL   <- shift.vector(Value5,-2, fill = NA)
  Value5L    <- shift.vector(Value5,-1, fill = NA)
  Value5R    <- shift.vector(Value5, 1, fill = NA)
  Value5RR   <- shift.vector(Value5, 2, fill = NA)
  
  # this is the funny KNN operation
  # B11 is central
  out  <-
    (-Value5RR + 4 * (Value5L + Value5R) + 10 * Value5 - Value5LL) / 16
  
  # cut back down (depending) and name
  names(out) <- Age5
  out
}

#' @export
#' @rdname smooth_age_5_un
united_nations_smth <- smooth_age_5_un

#' A strong method of population count smoothing
#' @description Smooth population counts in 5-year age groups.
#' @details The open age group is aggregated down to be evenly divisible by 10. This method accounts for the youngest and oldest age groups. Age classes must be cleanly groupable to 5-year age groups. All age classes are returned, but the strongest adjustment occurs in ages bounded by \code{ageMin} and \code{ageMax}. To be clear \code{ageMax} refers to the lower bound of the highest age class, inclusive. So, if you want a ceiling of 70 (default), specify 65. Counts are not constrained within this range, but the youngest 10-year age group and penultimate 10-year age group are perturbed but constrained to their original totals. The oldest 10-year age group is unaffected.
#' @param Value numeric vector of counts in single, abridged, or 5-year age groups.
#' @param Age numeric vector of ages corresponding to the lower integer bound of the counts.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}.
#' @param ageMin integer. The lowest age included included in intermediate adjustment. Default 10.
#' @param ageMax integer. The highest age class included in intermediate adjustment. Default 65.
#' @return numeric vector of smoothed counts in 5-year age groups.
#' @export
#' @examples
#' Ages         <- seq(0, 80, by = 5)
#' strongtest <- c(646617,511270,386889,317345,273736,240058,218645,188297,
#' 		153931, 124347,93254,71858,53594,39721,27887,18092,34729 )
#' strong_result <- smooth_age_5_strong(pop5m_pasex,Ages,TRUE)
#' # differences due to intermediate rounding in spreadsheet (bad practice IMO)
#' all(abs(strong_result - strongtest) < 1, na.rm = TRUE)
#' \dontrun{
#' plot(Ages, pop5m_pasex)
#' lines(as.integer(names(strong_result)),strong_result)
#' }
#' @references
#' \insertRef{carrier1959reduction}{DemoTools}
# plz add Demgen 1994 and PASEX. I can't download above citation while
# on vacation. plz check if this method is discussed there as well.

smooth_age_5_strong <- function(Value,
                        Age,
                        OAG = TRUE,
                        ageMin = 10,
                        ageMax = 65) {
  if (as.character(match.call()[[1]]) == "strong_smth") {
    warning("please use smooth_age_5_strong() instead of strong_smth().", call. = FALSE)
  }
  
  # these values are not used, it's just for lengths, and to make sure we
  # end on an even 10. Technically we could even provide data in 10-year
  # age groups and it'd still not break.
  Value1     <- graduate_uniform(Value = Value, Age = Age, OAG = OAG)
  Value5     <-
    groupAges(Value1, Age = as.integer(names(Value1)), N = 5)
  N          <- length(Value5)
  Age5       <- as.integer(names(Value5))
  
  # would need to move this up to ensure?
  # or in case of 85+ would we want to keep 80-84, 85+ as-is?
  Value10    <- groupAges(Value, Age = Age, N = 10)
  
  # what OAG is a strange digit? Then take OAG after grouping.
  if (OAG) {
    OAGvalue <- Value5[length(Value5)]
    Value10[length(Value10)] <- NA
    Value5[length(Value5)]   <- NA
  }
  Age5       <- as.integer(names(Value5))
  Age10      <- as.integer(names(Value10))
  
  # subtotal
  indsub     <- Age10 >= ageMin & Age10 <= ageMax
  SubTot     <- sum(Value10[indsub])
  
  #	# get staggered vectors
  Value10L    <- shift.vector(Value10,-1, fill = NA)
  Value10R    <- shift.vector(Value10, 1, fill = NA)
  # this is the funny KNN operation
  # B11 is central
  Value10Pert <- (Value10 * 2 + Value10L + Value10R) / 4
  Value10Pert[is.na(Value10Pert)] <- Value10[is.na(Value10Pert)]
  
  # rescale ages between min and max to sum to original
  Value10Adj  <- Value10Pert
  Value10Adj[indsub] <-
    Value10Adj[indsub] * SubTot / sum(Value10Adj[indsub])
  
  # again get staggered vectors
  V10adjLL    <- shift.vector(Value10Adj,-2, fill = NA)
  V10adjL     <- shift.vector(Value10Adj,-1, fill = NA)
  V10adjR     <- shift.vector(Value10Adj, 1, fill = NA)
  V10adjRR    <- shift.vector(Value10Adj, 2, fill = NA)
  
  # alternating calc, with differences at tails
  vevens     <- (-V10adjR + 11 * Value10Adj + 2 * V10adjL) / 24
  # tails different
  vevens[1]  <-
    (8 * Value10Adj[1] + 5 * V10adjL[1] - V10adjLL[1]) / 24
  lastind     <- which(is.na(vevens))[1]
  vevens[lastind] <-
    Value10Adj[lastind] - (8 * Value10Adj[lastind] + 5 * V10adjR[lastind] - V10adjRR[lastind]) / 24
  # odds are complement
  vodds      <- Value10Adj - vevens
  
  # prepare output
  # prepare output
  interleaf  <- c(rbind(vodds, vevens))
  # produce results vector
  out        <- Value5 * NA
  n          <- min(c(length(interleaf), N))
  out[1:n]   <- interleaf[1:n]
  
  # if OA ends in 5, then we can save penultimate value too
  na.i       <- is.na(out)
  out[na.i]  <- Value5[na.i]
  if (OAG) {
    out[N] <- OAGvalue
  }
  
  # what if OAis e.g. 85?
  out
}
#' @export
#' @rdname smooth_age_5_strong
strong_smth <- smooth_age_5_strong

#' G. Feeney's method of smoothing counts in 5-year age groups.
#' @description If age heaping is much worse on 0's than on 5's then even counts in 5-year age bins can preserve a sawtooth pattern. Most graduation techniques translate the zig-zag/sawtooth pattern to a wave pattern. It is not typically desired. This method redistributes counts 'from' every second 5-year age group in a specified range 'to' the adjacent age groups. How much to redistribute depends on a detection of roughness in the 5-year binned data, which follows the formulas recommended by Feeney. This method does not alter the total population count, counts in the youngest 10 ages, nor in old ages. 10-year age groups in the middle age range are not constrained.
#' @details This function calls \code{smooth_age_5_zigzag_inner()}, but prepares data in a way consistent with other methods called by \code{smooth_age_5()}. It is probably preferable to call \code{zigzag()} from the top level, or else call this method from \code{smooth_age_5()} for more control over tail imputations.
#' @param Value numeric vector of (presumably) counts in 5-year age groups.
#' @param Age integer vector of age group lower bounds.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}.
#' @param ageMin integer. Lower age bound to adjust values.
#' @param ageMax integer. Upper age bound to adjust values.
#' @return numeric vector of smoothed counts in 5-year age groups.
#' @export
#' @references
#' Feeney, G. 2013 "Removing "Zigzag" from Age Data," http://demographer.com/white-papers/2013-removing-zigzag-from-age-data/
#' @examples
#' Age <- c(0,1,seq(5,90,by=5))
#' # defaults
#' zz <- smooth_age_5_zigzag(dth5_zigzag, Age, OAG = TRUE, ageMin = 40, ageMax = 90)
#' \dontrun{
#' plot(Age, dth5_zigzag)
#' lines(as.integer(names(zz)),zz)
#' }
smooth_age_5_zigzag <- function(Value,
                        Age,
                        OAG = TRUE,
                        ageMin = 40,
                        ageMax = max(Age) - max(Age) %% 10 - 5) {
  
  if (as.character(match.call()[[1]]) == "zigzag_smth") {
    warning("please use smooth_age_5_zigzag() instead of zigzag_smth().", call. = FALSE)
  }
  # insist on 5-year age groups
  Value <- groupAges(Value, Age = Age, N = 5)
  Age   <- as.integer(names(Value))
  
  Smoothed <- smooth_age_5_zigzag_inner(
    Value = Value,
    Age = Age,
    OAG = OAG,
    ageMin = ageMin,
    ageMax = ageMax
  )
  # put NAs in for unadjusted ages,
  Smoothed[Smoothed == Value] <- NA
  
  Smoothed
}

#' @export
#' @rdname smooth_age_5_zigzag
zigzag_smth <- smooth_age_5_zigzag


#' Smooth in 5-year age groups using a moving average
#' @description Smooth data in 5-year age groups.
#' @details This function calls \code{smooth_age_5_zigzag_inner()}, but prepares data in a way consistent with other methods called by \code{smooth_age_5()}. It is probably preferable to call \code{zigzag()} from the top level, or else call this method from \code{agesmth()} for more control over tail imputations.
#' @param Value numeric vector of (presumably) counts in 5-year age groups.
#' @param Age integer vector of age group lower bounds.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}.
#' @param n integer. The width of the moving average. Default 3 intervals (x-5 to x+9).
#' @return numeric vector of smoothed counts in 5-year age groups.
#' @details This function calls \code{mav()}, which itself relies on the more general \code{ma()}. We lose the lowest and highest ages with this method, unless \code{n=1}, in which case data is returned in the original 5-year age groups. The total population count is not constrained to sum to the orignal total.
#' @examples
#' Age <- c(0,1,seq(5,90,by=5))
#' # defaults
#' ns   <- sapply(1:5,smooth_age_5_mav,Value=dth5_zigzag,Age=Age,OAG=TRUE)
#' cols <- paste0(gray(seq(.8,0,length=5)),"A0")
#' lwds <- seq(4,1,length=5)
#' \dontrun{
#' plot(Age, dth5_zigzag,pch=16)
#' matplot(as.integer(rownames(ns)),ns,type='l',
#' 		col = cols,
#' 		lty = 1,
#' 		add = TRUE,
#' 		lwd = lwds)
#' legend("topright", col = cols, lty = 1, lwd = lwds, legend = paste("n =",1:5))
#' }
#' @export

smooth_age_5_mav <- function(Value,
                     Age,
                     OAG = TRUE,
                     n = 3) {
  if (as.character(match.call()[[1]]) == "mav_smth") {
    warning("please use smooth_age_5_mav() instead of mav_smth().", call. = FALSE)
  }
  Value <- groupAges(Value, Age = Age, N = 5)
  Age   <- as.integer(names(Value))
  
  Smoothed <- mav(
    Value = Value,
    Age = Age,
    OAG = OAG,
    n = n
  )
  
  Smoothed
}

#' @export
#' @rdname smooth_age_5_mav
mav_smth <- smooth_age_5_mav


# Author: Juan Galeano
# handle OAG with care.
###############################################################################

#' Feeney'S formula on 9 years to correct for heaping on multiples of 5.

#' @description  Fenney's technique for correcting age distributions for heaping on multiples of five.

#' @param Value numeric. A vector of demographic counts in single age groups.
#' @param Age numeric or character. A vector with ages in single years.
#' @param maxit integer. Maximum number of iterations.
#' @param OAG logical. Is the final age group open? Default \code{FALSE}.

#' @details \code{Value} can be given in single or 5-year age groups.

#' @return a vector of adjusted counts in 5-year age groups
#'
#' @export
#' @references
#' \insertRef{feeney1979}{DemoTools}

#' @examples
#' # data from feeney1979, Table 1, page 12: Population of Indonesia, 22 provinces,
#' # by single year of age: Census of 24 September 1971.
#'  Pop        <- c(2337,3873,3882,3952,4056,3685,3687,3683,3611,3175,
#'          3457,2379,3023,2375,2316,2586,2014,2123,2584,1475,
#'          3006,1299,1236,1052,992,3550,1334,1314,1337,942,
#'          3951,1128,1108,727,610,3919,1221,868,979,637,
#'          3409,887,687,533,313,2488,677,426,524,333,
#'          2259,551,363,290,226,1153,379,217,223,152,
#'          1500,319,175,143,89,670,149,96,97,69,
#'          696,170,60,38,23,745)
#'  Ages       <- c(0:75)
#'  result     <- smooth_age_5_feeney(Pop, Ages, OAG = TRUE)
#'  A5         <- names2age(result)
#'  V5         <- groupAges(Pop,Ages)
#'  \dontrun{
#'  plot(Ages, Pop, type= 'l')
#'  segments(A5,
#'		  result/5,
#'		  A5+5,
#'		 result/5,
#' 		 col = "red")
#' segments(A5,
#'		 V5/5,
#'		 A5+5,
#'		 V5/5,
#'		 col = "blue")
#'  legend("topright",col=c("black","blue","red"),
#'    lty=c(1,1,1),
#'    legend=c("recorded 1","recorded 5","corrected 5"))
#' }

smooth_age_5_feeney          <- function(Value,
                                         Age,
                                         maxit = 200,
                                         OAG = FALSE) {
  if (as.character(match.call()[[1]]) == "T9R5L") {
    warning("please use smooth_age_5_feeney() instead of T9R5L().", call. = FALSE)
  }
  
  # ages need to be single to use this method.
  stopifnot(is_single(Age))
  TOT          <- sum(Value, na.rm = TRUE)
  
  # handle OAG with care
  if (OAG) {
    NN         <- length(Value)
    OAvalue    <- Value[NN]
    OA         <- Age[NN]
    Value      <- Value[-NN]
    Age        <- Age[-NN]
    
  }
  
  V5           <- groupAges(Value, Age, N = 5)
  A5           <- names2age(V5)
  i5           <- Age %% 5 == 0
  
  # TR: is this length-vulnerable?
  V50          <- Value[i5]
  V54          <- V5 - V50
  
  
  # need N anyway
  N            <- length(V5)
  # internal function used iteratively
  f_adjust     <- function(v5, v4) {
    N          <- length(v4)
    Bup        <- c(1, shift.vector(v4,-1, fill = 0) + v4)[1:N]
    FC         <- (8 / 9) * ((v5 + Bup) / Bup)
    FC[1]      <- 1
    POB1       <- c(v5 - (FC - 1) * Bup)
    POB2       <- (c(c(shift.vector(FC,-1, fill = 0) + FC)[-N],
                     (FC[N] + 1)) - 1) * v4
    # return a list of 3 components
    aj         <- list(FC, POB1, POB2)
    aj
  }
  
  
  # adjustment loop
  for (i in 1:maxit) {
    adjust     <- f_adjust(v5 = V50, v4 = V54)
    V50        <- adjust[[2]]
    V54        <- adjust[[3]]
    if (all(abs(adjust[[1]]) < 1e-8)) {
      break
    }
  }
  
  G            <- (V50 * .6)[c(2:(N - 1))]
  H            <- (V50 * .4)[c(3:N)]
  I            <- G + H #f(x+2.5)
  
  # corrected, but unknowns still need to be redistributed
  #I5          <- rescale_vector(I * 5,sum(V5[2:(N-1)]))
  I5           <- I * 5
  out          <- c(V5[1], I5, V5[N])
  
  if (OAG) {
    A5         <- c(A5, OA)
    out        <- c(out, OAvalue)
  }
  names(out)   <- A5
  # rescale to sum, inclusing open age group and boudned tails
  out          <- rescale_vector(out, TOT)
  out
}

#' @export
#' @rdname smooth_age_5_feeney
T9R5L <- smooth_age_5_feeney

#' Smooth populations in 5-year age groups using various methods
#'
#' @description Smooth population counts in 5-year age groups using the Carrier-Farrag,
#' Karup-King-Newton, Arriaga, United Nations, Stong, or Zigzag methods. Allows for imputation
#' of values in the youngest and oldest age groups for the Carrier-Farrag, Karup-King-Newton,
#' and United Nations methods.

#' @details The Carrier-Farrag, Karup-King-Newton (KKN), and Arriaga methods do not modify the totals
#' in each 10-year age group, whereas the United Nations, Strong, Zigzag, and moving average (MAV) methods do. The age intervals
#' of input data could be any integer structure (single, abridged, 5-year, other), but
#' output is always in 5-year age groups. All methods except the United Nations and MAV methods
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
#' compare the results from a variety of methods.
#'
#' @param Value numeric vector of counts in single, abridged, or 5-year age groups.
#' @param Age integer vector of ages corresponding to the lower integer bound of the counts.
#' @param method character string. Options include \code{"Carrier-Farrag"},\code{"Arriaga"},\code{"KKN"},\code{"United Nations"}, \code{"Strong"}, and \code{"Zigzag"}. See details. Default \code{"Carrier-Farrag"}.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}.
#' @param ageMin integer. The lowest age included included in intermediate adjustment. Default 10. Only relevant for Strong method.
#' @param ageMax integer. The highest age class included in intermediate adjustment. Default 65. Only relevant for Strong method.
#' @param n integer. The width of the moving average. Default 3 intervals (x-5 to x+9). Only relevant for moving average method.
#' @param young.tail \code{NA} or character. Method to use for ages 0-9. See details. Default \code{"original"}.
#' @param old.tail \code{NA} or character. Method to use for the final age groups. See details. Default \code{"original"}.
#'
#' @return numeric vector of smoothed counts in 5-year age groups.
#' @export
#'
#' @examples
#' Ages         <- seq(0, 80, by = 5)
#'
#' # names a bit flexible:
#' cf <- smooth_age_5(Value = pop5m_pasex,
#'		Age = Ages,
#'		method = "Carrier-Farrag",
#'		OAG = TRUE)
#'# old.tail be default same as young.tail
#'# "cf" also works
#'
#'# no need to specify tails for Arriaga or Strong
#'arr <- smooth_age_5(Value = pop5m_pasex,
#'		Age = Ages,
#'		method = "Arriaga",
#'		OAG = TRUE)
#'strong <- smooth_age_5(Value = pop5m_pasex,
#'		Age = Ages,
#'		method = "Strong",
#'		OAG = TRUE)
#'# other methods:
#'un <- smooth_age_5(Value = pop5m_pasex,
#'		Age = Ages,
#'		method = "United Nations",
#'		OAG = TRUE)
#'kkn <- smooth_age_5(Value = pop5m_pasex,
#'		Age = Ages,
#'		method = "KKN",
#'		OAG = TRUE)
#' # zigzag, not plotted.
#' zz <- smooth_age_5(pop5m_pasex,Ages,OAG=TRUE,method="Zigzag",ageMin = 30, ageMax = 80)
#' # mav, not plotted.
#' ma3 <- smooth_age_5(pop5m_pasex,Ages,OAG=TRUE,method="MAV",n=3)
#'
#'\dontrun{
#'	plot(Ages,pop5m_pasex,pch=16)
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
#'  Age <- 0:99
#'
#'  V5 <- groupAges(pop1m_pasex, Age=Age)
#'  Age5 <- as.integer(names(V5))
#'  cf2 <- smooth_age_5(Value = pop1m_pasex,
#'		  Age = Age,
#'		  method = "Carrier-Farrag",
#'		  OAG = TRUE)
#'  st2 <- smooth_age_5(Value = pop1m_pasex,
#'		  Age = Age,
#'		  method = "Strong",
#'		  OAG = TRUE)
#'  \dontrun{
#'  plot(Age,pop1m_pasex,pch=16)
#'  lines(Age,graduate_uniform(V5,Age=Age5,OAG=FALSE), lty=2, lwd = 2)
#'  lines(Age,graduate_uniform(cf2,Age=Age5,OAG=FALSE),col="blue")
#'  lines(Age,graduate_uniform(st2,Age=Age5,OAG=FALSE),col="red")
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
#'spr1 <- graduate_sprague(pop1m_pasex, Age=Age,OAG=FALSE)
#'spr2 <- graduate_sprague(cf2, Age=Age5,OAG=FALSE)
#'spr3 <- graduate_sprague(st2, Age=Age5,OAG=FALSE)
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
#'
smooth_age_5 <- function(Value,
                    Age,
                    method = c("Carrier-Farrag",
                               "KKN",
                               "Arriaga",
                               "United Nations",
                               "Strong",
                               "Zigzag",
                               "MAV"),
                    OAG = TRUE,
                    ageMin = 10,
                    ageMax = 65,
                    n = 3,
                    young.tail = c("Original", "Arriaga", "Strong", NA),
                    old.tail = young.tail) {
  if (as.character(match.call()[[1]]) == "agesmth") {
    warning("please use smooth_age_5() instead of agesmth().", call. = FALSE)
  }
  
  method <- match.arg(method, c("Carrier-Farrag",
                                "KKN",
                                "Arriaga",
                                "United Nations",
                                "Strong",
                                "Zigzag",
                                "MAV"))
  young.tail <- match.arg(young.tail, c("Original", "Arriaga", "Strong", NA))
  old.tail   <- match.arg(old.tail, c("Original", "Arriaga", "Strong", NA))
  
  method     <- simplify.text(method)
  young.tail <- simplify.text(young.tail)
  old.tail   <- simplify.text(old.tail)
 
  
   
  if (missing(Age)) {
    Age      <- as.integer(names(Value))
  }
  stopifnot(length(Age) == length(Value))
  # carrierfarrag or cf
  if (method %in% c("cf", "carrierfarrag")) {
    out      <- smooth_age_5_cf(Value = Value,
                               Age = Age,
                               OAG = OAG)
  }
  
  # stong
  if (method == "strong") {
    out <-
      smooth_age_5_strong(
        Value = Value,
        Age = Age,
        OAG = OAG,
        ageMin = ageMin,
        ageMax = ageMax
      )
  }
  
  # un or unitednations
  if (method %in% c("un", "unitednations")) {
    out <- smooth_age_5_un(Value = Value,
                               Age = Age,
                               OAG = OAG)
  }
  
  # arriaga
  if (method  == "arriaga") {
    out <- smooth_age_5_arriaga(Value = Value,
                        Age = Age,
                        OAG = OAG)
  }
  
  # kkn kkingnewton karupkingnewton
  if (method %in% c("kkn", "kkingnewton", "karupkingnewton")) {
    out <- smooth_age_5_kkn(Value = Value,
                    Age = Age,
                    OAG = OAG)
  }
  # TR: new Feeney method added July 31, 2018
  # TR: ruh roh, we have another feeney method T9R5L, now called feeney.
  if (method %in% c("zigzag")) {
    # however, need to make it so NAs returned in unaffected ages?
    # or make the user call it in various runs and graft together.
    out <- smooth_age_5_zigzag(
             Value = Value,
             Age = Age,
             OAG = OAG,
             ageMin = ageMin,
             ageMax = ageMax)
  }

  if (method %in% c("feeney")) {
    # however, need to make it so NAs returned in unaffected ages?
    # or make the user call it in various runs and graft together.
    out <- smooth_age_5_feeney(
             Value = Value,
             Age = Age,
             OAG = OAG,
             maxit = 200)
  }
  
  if (method %in% c("mav", "ma", "movingaverage")) {
    # however, need to make it so NAs returned in unaffected ages?
    # or make the user call it in various runs and graft together.
    out <- smooth_age_5_mav(
             Value = Value,
             Age = Age,
             OAG = OAG,
             n = n)
  }
  # -------------------------------
  # clean tails
  nas     <- is.na(out)
  if (any(nas) & (!is.na(old.tail) | !is.na(young.tail))) {
    nrle             <- rle(as.integer(nas))
    original         <- groupAges(Value, Age = Age, N = 5)
    arriaga          <- smooth_age_5_arriaga(Value, Age = Age, OAG = OAG)
    strong           <- smooth_age_5_strong(Value, Age = Age, OAG = OAG)
    # are the final entries NAs?
    if (nrle$values[length(nrle$values)] == 1 & !is.na(old.tail)) {
      nrle$values[1] <- 0
      old.ind        <- as.logical(rep(nrle$values, times = nrle$lengths))
      # do we want original values?
      if (old.tail %in% c("o", "orig", "original")) {
        stopifnot(length(original) == length(out))
        out[old.ind] <- original[old.ind]
      }
      # or arriaga?
      if (old.tail == "arriaga") {
        stopifnot(length(arriaga) == length(out))
        out[old.ind] <- arriaga[old.ind]
      }
      # or strong?
      if (old.tail == "strong") {
        stopifnot(length(strong) == length(out))
        out[old.ind] <- strong[old.ind]
      }
      
    }
    nrle             <- rle(as.integer(nas))
    # take care of young tail
    if (nrle$values[1] == 1 & !is.na(young.tail)) {
      nrle$values[length(nrle$values)] <- 0
      young.ind        <-
        as.logical(rep(nrle$values, times = nrle$lengths))
      
      if (young.tail %in% c("o", "orig", "original")) {
        stopifnot(length(original) == length(out))
        out[young.ind] <- original[young.ind]
      }
      # or arriaga?
      if (young.tail == "arriaga") {
        stopifnot(length(arriaga) == length(out))
        out[young.ind] <- arriaga[young.ind]
      }
      # or strong?
      if (young.tail == "strong") {
        stopifnot(length(strong) == length(out))
        out[young.ind] <- strong[young.ind]
      }
    }
  } # end tail operations
  
  out
}


#' @export
#' @rdname smooth_age_5
agesmth <- smooth_age_5
