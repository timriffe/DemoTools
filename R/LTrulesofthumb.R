
# Author: tim
###############################################################################

# introducing a couple rules of thumb based on HMD relationships and simple assumptions.
# these are intended to kick in only if a demographic process-informed method isn't
# feasible.

#' rule of thumb for splitting infants from deaths under 5
#' @description Given the log crude death rate for ages 0-5, \eqn{log(M(0,5)} there is a rough 2-slope linear relationship with the log of the proportion of deaths under 5 that occur in the first year of life. This rule uses this two-slope linear relationship to estimate the proportion of deaths in age 0, and then split these from deaths under 5. The relationship is slightly different between males and females. This method is only to be invoked if neither population nor deaths are available in a separate tabulation for single age 0.
#' @details This is an elsewhere-undocumented relationship derived from the whole of the HMD. We used the \code{segmented} package to fit a 2-slope linear model. This can (and should) be reproduced using data from a more diverse collection, and even as-is the data should be subset to only those observations where deaths and populations were not split using HMD methods. You can reproduce the analysis given a data set in the format shown (but not executed) in the examples and following the code steps shown there.
#' 
#' Regarding argument specification, \code{D04} is required, in which case either \code{M04} or \code{P04} can be given to continue.
#' @param D04 numeric. Deaths under age 5.
#' @param M04 numeric. Death rate under age 5.
#' @param P04 numeric. Exposure under age 5. A population estimate will do.
#' @param Sex character, either \code{"m"} or \code{"f"}.
#' @return Estimated deaths in age 0
#' @export
#' @references 
#' \insertRef{Vito2008}{DemoTools}
#' 
#' Human Mortality Database
#' @examples
#' 
#' # to reproduce the coefficient estimation 
#' # that the method is based on:
#' \dontrun{
#' # get data in this format:
#' 	# dput(head(Dat))
#' Dat <- structure(list(
#' 		lDR0 = c(-0.182459515740971, -0.147321312595521, 
#' 		         -0.138935222455608, -0.156873832361186,
#' 				 -0.134782094278661, -0.135213007819874), 
#' 		lM5 = c(-4.38578221044468, -4.56777854985643, 
#' 				-4.58851248767896, -4.57684656734645, 
#' 				-4.62854681692805, -4.61294106314254)), 
#'         .Names = c("lDR0", "lM5"), 
#'         class = c("data.frame"), 
#' 		row.names = c(NA, -6L))
#'  # where lDRO is log(D0 / D0_4)
#'  # i.e. lof of proportion of deaths < 5 in first year of life
#'  # and lM5 is log(M0_4) 
#'  # i.e. log of death rate in first 5 years of life
#' 
#'  # then first fit a linear model:
#' 	obj  <- lm(lDR0~lM5, data = Dat)
#'  # use segmented package:
#' 	seg  <- segmented::segmented(obj)
#'  # breakpoint:
#' 	seg$psi[2]     # brk
#'  # first intercept:
#' 	seg$coef[1]    # int1
#'  # first slope:
#' 	seg$coef[2]    # s1
#'  # difference in slope from 1st to second:
#' 	seg$coef[3]    # ds1
#'  # make Dat come from some other dataset and you'll get different coefs,
#'  # it'd be possible to have these in families maybe, and in any case 
#'  # different for males and females. This is just a rough start, to be 
#'  # replaced if someone offers a superior method. These 
#' 
#' }
#' D0_4 <- 2e4
#' M0_4 <- 5/1000
#' P0_4 <- 4e6
#' # function usage straightforward, also vectorized.
#' D0   <- M04_2_D0(D0_4, M0_4, Sex = "m")
#' # deaths in ages 1-4 are a separate step.
#' D1_4 <- D0_4 - D0
#' 
#' # to get M0_4 it's best to follow these steps:
#' # 1) get M0 using M04_2_M0(M0_4),
#' M0   <- M04_2_M0(M0_4)
#' # 2) get denom using P0 = D0 / M0
#' P0   <- D0 / M0
#' # 3) get denom P1_4  as P0_4 - P0
#' P1_4 <- P0_4 - P0
#' # 4) M1_4 = D1_4 / P1_4.
#' M1_4 <- D1_4 / P1_4
#' 
#' \dontrun{
#' plot(NULL, type = 'n', xlim = c(0, 5), ylim = c(1e-3, .025), log = "y",
#' 		xlab = "Age", ylab = "log(rate)")
#' segments(0, M0_4, 5, M0_4)
#' segments(0, M0, 1, M0)
#' segments(1, M1_4, 5, M1_4)
#' text(1, c(M0, M1_4, M0_4), c("M0", "M1_4", "M0_4"), pos = 3)
#' }
M04_2_D0 <- function(D04, M04, P04, Sex = c("m","f")[1]){
	if (missing(M04)){
		M5 <- D04 / P04
	} else {
		M5 <- M04
	}
	
	lM5   <- log(M5)
	
	if (length(Sex)  == 1 & length(lM5) > 1){
		Sex   <- rep(Sex,length(lM5))
	}
	# get coefs, from segmented regression: 
	coefs <- structure(c(-4.69468169447114, -0.101269561161854, 0.0171297174685972, 
					-0.194795399634278, -4.38840457546843, -0.0800885256382603, 0.0208247842765569, 
					-0.192054558152916), .Dim = c(4L, 2L), .Dimnames = list(c("brk", 
							"int1", "s1", "ds1"), c("f", "m")))
	s2     <- coefs["s1",] + coefs["ds1",]
	int2   <- coefs["int1",] + coefs["s1",] * coefs["brk",] - s2 * coefs["brk",]
	coefs  <- rbind(coefs, s2, int2)
	# based on a segmented regression of log(1D0/5D0)~log(5M0)
	lrat  <- ifelse(Sex == "m",
			# if male
			ifelse(lM5 < coefs["brk","m"],
					coefs["int1","m"] + lM5 * coefs["s1","m"],
					coefs["int2","m"] + lM5 * coefs["s2","m"]),
			# if female
			ifelse(lM5 < coefs["brk","f"],
					coefs["int1","f"] + lM5 * coefs["s1","f"],
					coefs["int2","f"] + lM5 * coefs["s2","f"])
	)
	# back to ratio
	rat   <- exp(lrat)
	prop0 <- 1 / (1 + rat)
	# return D0
	prop0 * D04
}

#' rule of thumb for estimating infant mortality rate from under 5 mortality
#' @description Given the log crude death rate for ages 0-5, \eqn{log(M(0,5)} there is a strong 2-slope linear relationship with the log of the infant death rate. With this method th euser supplies the under 5 death rate and the infant death rate is returned. This method could be invoked if either population or (exclusive or) deaths are tabulated for infants, or if only the under 5 mortality rate is available
#' @details This is an elsewhere-undocumented relationship derived from the whole of the HMD. We used the \code{segmented} package to fit a 2-slope linear model. This can (and should) be reproduced using data from a more diverse collection, and even as-is the data should be subset to only those observations where deaths and populations were not split using HMD methods. You can reproduce the analysis given a data set in the format shown (but not executed) in the examples and following the code steps shown there.
#' 
#' Regarding argument specification, either \code{M04} *or* \code{D04} and \code{P04} can both be given.
#' @param D04 numeric. Deaths under age 5.
#' @param M04 numeric. Death rate under age 5.
#' @param P04 numeric. Exposure under age 5. A population estimate will do.
#' @param Sex character, either \code{"m"} or \code{"f"}.
#' @return Estimated deaths in age 0
#' @export
#' @references 
#' \insertRef{Vito2008}{DemoTools}
#' 
#' Human Mortality Database
#' @examples
#' 
#' # to reproduce the coefficient estimation 
#' # that the method is based on:
#' \dontrun{
#' # get data in this format:
#' 	# dput(head(Dat))
#' Dat <- structure(
#'		list(
#'		  lm0 = c(-2.92192434640927, -3.06165842367016, 
#'				  -3.10778177736261, -3.14075804560425, 
#'				  -3.20005407062761, -3.22489389719376
#'				), 
#'		  lM5 = c(-4.38578221044468, -4.56777854985643, 
#'				  -4.58851248767896, -4.57684656734645, 
#'				  -4.62854681692805, -4.61294106314254)), 
#'        .Names = c("lm0", "lM5"), 
#'		class = c("data.frame"), 
#'		row.names = c(NA, -6L))
#'  # where lm0 is log(M0)
#'  # i.e. log of infant death rate
#'  # and lM5 is log(M0_4) 
#'  # i.e. log of death rate in first 5 years of life
#' 
#'  # then first fit a linear model:
#' 	obj  <- lm(lm0~lM5,data=Dat)
#'  # use segmented package:
#' 	seg  <- segmented::segmented(obj)
#'  # breakpoint:
#' 	seg$psi[2]     # brk
#'  # first intercept:
#' 	seg$coef[1]    # int1
#'  # first slope:
#' 	seg$coef[2]    # s1
#'  # difference in slope from 1st to second:
#' 	seg$coef[3]    # ds1
#'  # make Dat come from some other dataset and you'll get different coefs,
#'  # it'd be possible to have these in families maybe, and in any case 
#'  # different for males and females. This is just a rough start, to be 
#'  # replaced if someone offers a superior method. These 
#' 
#' }
#'
#' M0_4 <- 5/1000
#' M0   <- M04_2_M0(M0_4)
#' 
#' # M0 from this relationship is more reliable than other methods
#' # of independently splitting D0 or P0. So, if you're going to be
#' # splitting counts, then make use of M0 to force numerators and 
#' # denominators to conform to this estimate.
#' D0_4 <- 2e4
#' P0_4 <- 4e6
#' # function usage straightforward, also vectorized.
#' D0   <- M04_2_D0(D0_4, M0_4, Sex = "m")
#' # deaths in ages 1-4 are a separate step.
#' P0   <- D0 / M0
#' P1_4 <- P0_4 - P0
#' D1_4 <- D0_4 - D0
#' M1_4 <- D1_4 / P1_4
#' # and now we have all the pieces such that rate
#' # estimates conform.
#' 
#' \dontrun{
#' plot(NULL, type = 'n', xlim = c(0, 5), ylim = c(1e-3, .025), log = "y",
#' 		xlab = "Age", ylab = "log(rate)")
#' segments(0, M0_4, 5, M0_4)
#' segments(0, M0, 1, M0)
#' segments(1, M1_4, 5, M1_4)
#' text(1, c(M0, M1_4, M0_4), c("M0", "M1_4", "M0_4"), pos = 3)
#' }

M04_2_M0 <- function(M04, D04, P04, Sex = c("m","f")[1]){
	if (missing(M04)){
		M5 <- D04 / P04
	} else {
		M5 <- M04
	}
	
	lM5   <- log(M5)
	
	# handle lengths
	if (length(Sex)  == 1 & length(lM5) > 1){
		Sex   <- rep(Sex,length(lM5))
	}
	
	# ------------------------------
	# get coefs created
	coefs  <- structure(c(-4.51388967799889, 1.49906532702989, 1.01391932870931, 
					-0.226555050620318, -4.76819124030249, 1.47765861943473, 1.01019112250369, 
					-0.227872879237884), .Dim = c(4L, 2L), .Dimnames = list(c("brk", 
							"int1", "s1", "ds1"), c("m", "f")))
	s2     <- coefs["s1",] + coefs["ds1",]
	int2   <- coefs["int1",] + coefs["s1",] * coefs["brk",] - s2 * coefs["brk",]
	coefs  <- rbind(coefs, s2, int2)
	# ------------------------------
	
	# based on a segmented regression of log(1D0/5D0)~log(5M0)
	lm0  <- ifelse(Sex == "m",
			# if male
			ifelse(lM5 < coefs["brk","m"],
					coefs["int1","m"] + lM5 * coefs["s1","m"],
					coefs["int2","m"] + lM5 * coefs["s2","m"]),
			# if female
			ifelse(lM5 < coefs["brk","f"],
					coefs["int1","f"] + lM5 * coefs["s1","f"],
					coefs["int2","f"] + lM5 * coefs["s2","f"])
	)
	
	# back to rate
	exp(lm0)
}