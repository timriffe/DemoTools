
# Author: tim
###############################################################################
#


# OAG: logical is the last value an open age group that shoud be preserved?
# ... optional arguments to pass to \code{loess()}. notably \code{span} might be interesting, as it controls smoothness
#' wrapper to loess
#' @description smooth data over age, preserving the open age group if necessary. This is a simple wrapper to \code{stats::loess()} but using standard demographic arguments.
#' @details The total sum of \code{Value} is preserved in the output. One can control smoothness using the \code{spar} argument of \code{stats::loess()}. See \code{\link[stats]{loess}} for more details.
#' @inheritParams agesmth
#' @param ... optional arguments passed to \code{stats::loess()}
#' @export
#' @importFrom stats loess
#' @examples 
#' \dontrun{
#' Age <- 0:99
#' plot(Age,pop1m_pasex)
#' lines(Age, ogive_loess(pop1m_pasex, Age, FALSE))
#' }
#' @seealso \code{\link[stats]{loess}}
ogive_loess <- function(Value, Age, OAG = TRUE, ...){
	scale      <- sum(Value, na.rm = TRUE)
	N          <- length(Value)
	stopifnot(N == length(Age))
	
	age        <- Age
	
	# separate Open age group if desired.
	if (OAG){
		OA     <- Value[N]
		Value  <- Value[-N]
		age    <- age[-N]
		scale  <- scale - OA
	}
	
	fit        <- loess(Value ~ age, ...)
	out        <- fit$fitted
	
	# negatives not allowed!
	out[out < 0] <- 0
	# enforce sums
	out        <- rescale.vector(out,  scale) 
	
	# tack open age group back on
	if (OAG){
		out    <- c(out, OA)
	}
	
	# name vector
	names(out) <- Age
	
	out
}


# ---------
# degree integer degree of polynomial, default 2
# trans character. tranformation to Value prior to fitting? \code{"log"} is the
# only implemented possibility at this time. Leave empty otherwise.
# OAG: logical is the last value an open age group that shoud be preserved?
# ... other optional arguments to pass to \code{lm()}

#' fit a polynomial to demographic data
#' @description smooth data over age, preserving the open age group if necessary. This is a wrapper to \code{stats::lm()} but using standard demographic arguments.
#' @details The total sum of \code{Value} is preserved in the output. One can control smoothness by adjusting the degree of the polynomial (higher degree more wiggly). See \code{\link[stats]{lm}} for more details. One may wish to log transform the data before fitting the polynomial, in which case \code{trans} should be specified as \code{"log"}. \code{"power"} is also an option in which case the data is transformed by \code{Value^(1/pow)} before fitting, and then the prediction is back-transformed (negative values not allowed)-- This is friendlier in the case of 0s. For the log transformation, 0s have no weight.
#' @inheritParams agesmth
#' @param degree integer degree of polynomial. Default 2.
#' @param trans if a transformation is desired, either specify \code{"log"} or \code{"power"}, otherwise leave missing. 
#' @param pow if \code{"power"} specified for \code{trans} then the power transformation to apply.
#' @param ... optional arguments passed to \code{stats::lm()}
#' @export
#' @examples 
#' \dontrun{
#' Age <- 0:99
#' plot(Age,pop1m_pasex)
#' lines(Age, ogive_poly(pop1m_pasex, Age, OAG = FALSE))
#' lines(Age, ogive_poly(pop1m_pasex, Age, degree = 3, OAG = FALSE), col = "red")
#' 
#' }
#' @seealso \code{\link[stats]{lm}}
#' @importFrom stats lm


ogive_poly <- function(Value, Age, degree = 2, trans, pow = 2, OAG = TRUE, ...){
	scale     <- sum(Value, na.rm = TRUE)
	N         <- length(Value)
	stopifnot(N == length(Age))
	
	# rename because we modify
	age       <- Age
	value     <- Value
	
	# default regression weights
	w         <- rep(1, length(value))
	# separate Open age group if desired.
	if (OAG){
		OA    <- value[N]
		value <- value[-N]
		age   <- age[-N]
		scale <- scale - OA
		w     <- w[-N]
	}
	
	# -------------------------
	# log transform to constrain positive, if desired
	# indicator
	lTF       <- FALSE
	pTF       <- FALSE
	if (!missing(trans)){
		if (trans == "log"){
			lTF   <- TRUE
			value <- log(value)
			w[is.infinite(value)] <- NA
		} 
		if (trans == "power"){
			pTF <- TRUE
			value <- value ^ (1 / pow)
		}	
	}
	
	# -------------------------
	# build up polynomial expression
	polys     <- 2:degree
	polys2    <- paste0("age^", polys)
	polys3    <- paste0(" + I(",polys2,")")
	# final formula
	expr      <- paste0("Value. ~ age", paste(polys3, collapse = ""))
	
	# stick elements into ad hoc data.frame
	dataf     <- data.frame(Value. = value, age = age, w = w)
	# fit linear model
	fit       <- lm(expr, weight = w, data = dataf, ...)
	# get predicted values
	out       <- predict(fit)
	
	# -------------------------
	# transform back, if necessary
	if (lTF){
		out <- exp(out)
	} 		
	if (pTF){
		out <- out ^ pow
	}
	# -------------------------
	# closing steps
    # negatives not allowed!
    out[out < 0] <- 0
	# enforce sums
	out        <- rescale.vector(out,  scale) 
	
	# tack OAG back on, if necessary
	if (OAG){
		out    <- c(out, OA)
	}
	
	# name vector
	names(out) <- Age
	
	out
}


# same args as above:

ogive <- function(Value, Age, method = "loess", OAG = TRUE, degree = 2, trans, ...){
	
	if (method == "loess"){
		out <- ogive_loess(Value = Value, Age = Age, OAG = OAG, ...)
	}
	if (method == "poly"){
		out <- ogive_poly(Value = Value, Age = Age, degree = degree, trans = trans, OAG = OAG, ...)
	}
	out
}


