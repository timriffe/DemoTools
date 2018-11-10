# Author: sean
# edited by TR 9-Dec-2017 & 28 Feb-2018
###############################################################################

#' Interpolate between two population age distributions.
#' @description The interpolation is done by age (not cohort) using a linear or exponential function. This comes from the PAS spreadsheet called ADJINT.

#' @param Pop1   numeric. A vector of demographic counts by age at time 1 (earlier date).
#' @param Pop2   numeric. A vector of demographic counts by age at time 2 (later date).
#' @param Date1 date. The date corresponding to the population age distribution at time 1. See details for ways to express it.
#' @param Date2 date. The date corresponding to the population age distribution at time 2. See details for ways to express it.
#' @param DesiredDate date. The desired date of the output population age distribution. See details for ways to express it.
#' @param method string. The method to use for the interpolation, either "linear" or "exponential". Default "linear".
#' @param roundoutput logical. Whether or not to return integers. Default \code{FALSE}. 
#' @details The age group structure of the output is the same as that of the input. Ideally, \code{DesiredDate} should be between the Date1 and Date2. 
#' Dates can be given in three ways 1) a \code{Date} class object, 2) an unambiguous character string in the format \code{"YYYY-MM-DD"}, 
#' or 3) as a decimal date consisting in the year plus the fraction of the year passed as of the given date.
#' @return A vector of the interpolated population for the requested date.
#' @author Sean Fennel
#' @references 
#' \insertRef{PAS}{DemoTools}
#' @export
#' 
#' @examples 
#' # YYYY-MM-DD dates as character
#' EarlyDate       <- "1980-04-07"
#' LaterDate       <- "1990-10-06"
#' DesiredDate     <- "1985-07-01"
#' 
#' interpolatePop(popA_earlier, popA_later, EarlyDate, LaterDate, DesiredDate)
#' interpolatePop(popA_earlier, popA_later, EarlyDate, LaterDate, DesiredDate, method = "exponential")
#' \dontrun{
#'plot(popA_earlier, t ='l', col='red',
#'     ylim = c(400,1220000), 
#'     ylab = 'The counts', xlab = 'Age groups')
#'
#'lines(interpolatePop(popA_earlier, popA_later, EarlyDate, LaterDate, DesiredDate),
#'      t = 'l', col='black')
#'
#'lines(popA_later,
#'      t = 'l', col='blue')
#'
#'legend(12,1000000, 
#'       legend = c('Pop in 1980-04-07',
#'                  'Pop in 1985-07-01',
#'                  'Pop in 1990-10-06'),
#'       col=c('red','black','blue'),
#'       lty = 1)
#' }

interpolatePop <- function(Pop1, Pop2, Date1, Date2, DesiredDate, method = "linear", roundoutput = FALSE){
  
  stopifnot(length(Pop1) == length(Pop2))
	
  earlyDateDec      <- dec.date(Date1)
  laterDateDec      <- dec.date(Date2)
  desireDateDec     <- dec.date(DesiredDate)
  
  interpolateFactor <- (desireDateDec - earlyDateDec) / (laterDateDec - earlyDateDec)
  
  if (method == "exponential"){
    adjustedPop   <- Pop1 * exp(interpolateFactor * log(Pop2 / Pop1))
  }
  else {
    adjustedPop   <- Pop1 + (Pop2 - Pop1) * interpolateFactor
  }
  if (roundoutput){
	  adjustedPop <- round(adjustedPop)
  }
  return(adjustedPop)
}

# TR: 2018-08-20 either complementing or replacing previous pop interpolater.
# this one more general.
###############################################################################

#' Interpolate between two population age distributions.
#' @description The interpolation is done by age using a linear, exponential, or power function. This comes from the PAS spreadsheet called \code{AGEINT}. Be aware that this is not cohort-component interpolation.

#' @param popmat numeric. An age-period matrix of data to interpolate over. Age in rows, time in columns.
#' @param datesIn vector of dates. The exact dates in each column. See details for ways to express it.
#' @param datesOut vector of dates. The desired dates to interpolate to. See details for ways to express it.
#' @param method string. The method to use for the interpolation, either \code{"linear"}, \code{"exponential"}, or \code{"power"}. Default \code{"linear"}.
#' @param power numeric power to interpolate by, if \code{method = "power"}. Default 2.
#' @param ... arguments passed to \code{stats::approx}. For example, \code{rule}, which controls extrapolation behavior.
#' @details The age group structure of the output is the same as that of the input. Ideally, \code{datesOut} should be within the range of \code{datesIn}. If not, the left-side and right-side output are held constant outside the range if \code{rule = 2} is passed in, otherwise \code{NA} is returned (see examples). Dates can be given in three ways 1) a \code{Date} class object, 2) an unambiguous character string in the format \code{"YYYY-MM-DD"}, or 3) as a decimal date consisting in the year plus the fraction of the year passed as of the given date. 
#' 
#' For methods \code{"exponential"} and \code{"power"}, calculations are carried out using linear interpolation through \code{log(popmat)}, or \code{popmat^(1/power)} respectively, then back-transformed. If the data contain 0s, \code{method = "exponential"} will fail, but \code{"power"} would still generally work. 
#' 
#' @return numeric matrix (age-period) (or vector if \code{length(datesOut) == 1} of the interpolated data for the requested dates. 
#' @author Tim Riffe
#' @references 
#' \insertRef{PAS}{DemoTools}
#' @export
#' @importFrom stats approx
#' 
#' @examples
#' # Example of interpolating over time series of age-structured
#' # data. NOTE: age classes must conform. They don't have to be
#' # uniform, just the same in each year.
#' popmat <- structure(c(2111460L, 1971927L, 1651421L, 1114149L, 921120L, 
#' 859806L, 806857L, 847447L, 660666L, 567163L, 493799L, 322936L, 
#' 320796L, 163876L, 133531L, 120866L, 2548265L, 2421813L, 2581979L, 
#' 2141854L, 1522279L, 1321665L, 1036480L, 1024782L, 935787L, 789521L, 
#' 719185L, 481997L, 479943L, 268777L, 202422L, 167962L, 3753848L, 
#' 3270658L, 2930638L, 2692898L, 2222672L, 1788443L, 1514610L, 1491751L, 
#' 1054937L, 972484L, 796138L, 673137L, 554010L, 352264L, 293308L, 
#' 195307L, 3511776L, 3939121L, 4076601L, 3602857L, 2642620L, 2105063L, 
#' 1993212L, 1914367L, 1616449L, 1408498L, 994936L, 776537L, 706189L, 
#' 507085L, 315767L, 240302L, 3956664L, 3936424L, 3995805L, 4371774L, 
#' 4012982L, 3144046L, 2406725L, 2301514L, 2061252L, 1871502L, 1538713L, 
#' 1210822L, 896041L, 638954L, 401297L, 375648L), .Dim = c(16L, 
#' 5L), .Dimnames = list(c("0", "5", "10", "15", "20", "25", "30", 
#' "35", "40", "45", "50", "55", "60", "65", "70", "75"), c("1960", 
#' "1976", "1986", "1996", "2006")))
#' 
#' # We have irregularly spaced census inputs.
#' \dontrun{
#' 	matplot(seq(0,75,by=5),popmat,type='l',col=1:5)
#' 	text(20,popmat["20",],colnames(popmat),col=1:5)
#' }
#' # census dates as decimal: no need for year rounding.
#' # could also specify as "YYYY-MM-DD"
#' # (rounded to 2 digits here- not necessary)
#' datesIn  <- c(1960.73, 1976.90, 1986.89, 1996.89, 2006.89)
#' # could ask for single years or anything
#' datesOut <- seq(1960, 2005, by = 5)
#' 
#' # NOTE: rule is an argument to approx(), which does the work.
#' # if not passed in, 1960 is returned as NAs (outside the range
#' # of dates given). In this case, rule = 2 assumes constant 
#' # equal to nearest neighbor, meaning 1960 data will be returned
#' # as equal to 1960.73 input data. So, function should only be 
#' # allowed to extrapolate for short distances.
#' interp(popmat, datesIn, datesOut, "linear", rule = 2)
#' interp(popmat, datesIn, datesOut, "exponential", rule = 2)
#' interp(popmat, datesIn, datesOut, "power", rule = 2)
#' 
#' # ----------------------------------------------------------
#' # standard example data, taken from AGEINT PAS spreadsheet.
#' 
#'  # YYYY-MM-DD dates as character
#'  d1              <- "1980-04-07"
#'  d2              <- "1990-10-06"
#'  dout            <- "1985-07-01"
#'  
#'  (p_lin <- interp(cbind(popA_earlier,popA_later),c(d1,d2),dout,method = "linear"))
#'  (p_exp <- interp(cbind(popA_earlier,popA_later),c(d1,d2),dout,method = "exponential"))
#'  (p_pow <- interp(cbind(popA_earlier,popA_later),c(d1,d2),dout,method = "power"))
#'  
#'  \dontshow{
#'  print(dec.date(d1),digits=18)
#'  print(dec.date(d2),digits=18)
#'  print(dec.date(dout),digits=18)
#'  # these were set to different dates: recreate above dates...
#' expcheck <- c(142637.629438722,658772.565488024,881802.276313982,
#' 790568.050982602,631167.059816327,523673.857092558,426483.092945892,
#' 352375.329549735,349640.068841387,315083.384096823,243376.038026844,
#' 210990.890321914,179772.749510111,149358.515383768,112481.945265698,
#' 75813.0628150498,43827.7491909366,48878.656748965)
#' lincheck <- c(151295.597805253,698759.433295475,935326.835323442,
#' 838554.779337049,669478.29970015,555460.361866426,452370.21082717,
#' 373764.176717154,370862.888286807,334208.645663947,258148.731927464,
#' 223797.836576011,190684.784325522,158424.435131266,119309.492300436,
#' 80414.843580795,46488.0518583646,51845.5446984939)
#' 
#' stopifnot(all(abs(expcheck - p_exp)<1e-6))
#' stopifnot(all(abs(lincheck - p_lin)<1e-6))
#' # do controlled spreadsheet test, passing in decimal dates, etc
#' }
#' \dontrun{
#' age <- c(0,1,seq(5,80,by=5))
#' plot(age, popA_later, type = 'l', lwd = 2, main = "Interpolation options")
#' lines(age, popA_earlier, lwd = 2, lty = "8282")
#' lines(age, p_lin,col = "blue")
#' lines(age, p_pow,col = gray(.4))
#' lines(age, p_exp,col = "red")
#' text(30, popA_earlier[8], d1, pos = 1, cex = 1)
#' text(20, popA_later[6], d2, pos = 4, cex = 1)
#' legend("topright", lty = 1, col = c("blue", gray(.4), "red"),
#' 		legend = c("linear", "power 2", "exponential"), 
#' 		title = paste0("Interpolated 1985-07-01"))
#' 
#' }

interp <- function(popmat, 
		datesIn, 
		datesOut, 
		method = c("linear","exponential","power")[1], 
		power = 2,
		...){ # ... args passed to stats::approx . Can give control over extrap assumptions
	# a basic check
	stopifnot(ncol(popmat) == length(datesIn))
	
	# no sense documenting this wrapper ...
	.approxwrap <- function(x, y, xout, ...){
		stats::approx(x = x, y = y, xout = xout, ...)$y
	}
	
	# -----------------------
	# clean method declaration
	method <- tolower(method)
	if (grepl(method, pattern = "exp")){
		pattern <- "exponential"
	}
	if (grepl(method, pattern = "lin")){
		pattern <- "linear"
	}
	if (grepl(method, pattern = "pow")){
		pattern <- "power"
	}
	# -----------------------
	
	# coerce dates to decimal if necessary
	datesIn  <- sapply(datesIn, dec.date)
	datesOut <- sapply(datesOut, dec.date)
	
	
	# carry out transform 1
	if (method == "exponential"){
		if (any(popmat == 0)){
			stop("popmat contains 0s, so exponential interpolation won't work.\n
							Either handle 0s then do this, or\n
							maybe try method = 'power'.")
		}
		popmat <- log(popmat)
	}
	if (method == "power"){
		popmat <- popmat ^ (1 / power)
	}
	
	int <- apply(
			popmat,
			1,
			.approxwrap,
			x = datesIn,
			xout = datesOut,
			...)
	dims <- dim(int)
	if (!is.null(dims)){
		int <- t(int)
		rownames(int) <- rownames(popmat)
		colnames(int) <- datesOut
	} else {
		names(int)    <- rownames(popmat)
	}
	
	# transform back
	if (method == "exponential"){
		int <- exp(int)
	}
	if (method == "power"){
		int <- int^power
	}
	
	int
}

