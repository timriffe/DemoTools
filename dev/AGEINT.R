# Author: sean
# edited by TR 9-Dec-2017 & 28 Feb-2018
# TODO: make wrapper to fill in interpolations with various input years and output years.
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

#pop5     <- read.csv("/home/tim/Desktop/pop5.csv")
#pop5     <- pop5[,c("AgeStart","AgeSpan", "TimeMid","SexName", "DataValue")]
#pop5     <- pop5[pop5$AgeStart >= 0 & pop5$SexName == "Male", ]
#
#pop5     <- pop5[with(pop5, order(TimeMid,AgeStart)),]
#pop5     <- pop5[!(pop5$AgeStart == 0 & pop5$AgeSpan == -1) & pop5$AgeStart <= 75, ]
#pop52006 <- pop5[pop5$TimeMid > 2006,]
#p1_4     <- sum(pop52006[pop52006$AgeSpan == 1 & pop52006$AgeStart < 5 & pop52006$AgeStart > 0,"DataValue"])
#pop52006[pop52006$AgeSpan == 1 &pop52006$AgeStart == 1, "DataValue"] <- p1_4
#pop52006[pop52006$AgeSpan == 1 &pop52006$AgeStart == 1, "AgeSpan"] <- 4
#pop52006 <- pop52006[!(pop52006$AgeStart > 0 & pop52006$AgeSpan == 1),]
#pop52006 <- pop52006[-1, ]
#
#pop5 <- rbind(pop5[pop5$TimeMid < 2006, ], pop52006)
#
#
#ind1 <- pop5$AgeStart == 0 & pop5$AgeSpan == 1
#ind2 <- pop5$AgeStart == 1 & pop5$AgeSpan == 4
#pop5[ind1, "DataValue"] <- pop5[ind1, "DataValue"] +pop5[ind2, "DataValue"] 
#pop5[ind1, "AgeSpan"] <- 5
#pop5 <- pop5[!ind2, ]
#library(reshape2)

# Example of interpolating over time series of age-structured
# data. NOTE: age classes must conform. They don't have to be
# uniform, just the same in each year.
popmat <- structure(c(2111460L, 1971927L, 1651421L, 1114149L, 921120L, 
				859806L, 806857L, 847447L, 660666L, 567163L, 493799L, 322936L, 
				320796L, 163876L, 133531L, 120866L, 2548265L, 2421813L, 2581979L, 
				2141854L, 1522279L, 1321665L, 1036480L, 1024782L, 935787L, 789521L, 
				719185L, 481997L, 479943L, 268777L, 202422L, 167962L, 3753848L, 
				3270658L, 2930638L, 2692898L, 2222672L, 1788443L, 1514610L, 1491751L, 
				1054937L, 972484L, 796138L, 673137L, 554010L, 352264L, 293308L, 
				195307L, 3511776L, 3939121L, 4076601L, 3602857L, 2642620L, 2105063L, 
				1993212L, 1914367L, 1616449L, 1408498L, 994936L, 776537L, 706189L, 
				507085L, 315767L, 240302L, 3956664L, 3936424L, 3995805L, 4371774L, 
				4012982L, 3144046L, 2406725L, 2301514L, 2061252L, 1871502L, 1538713L, 
				1210822L, 896041L, 638954L, 401297L, 375648L), .Dim = c(16L, 
				5L), .Dimnames = list(c("0", "5", "10", "15", "20", "25", "30", 
						"35", "40", "45", "50", "55", "60", "65", "70", "75"), c("1960", 
						"1976", "1986", "1996", "2006")))
# census dates as decimal 
# (rounded to 2 digits here- not necessary)
datesIn  <- c(1960.73, 1976.90, 1986.89, 1996.89, 2006.89)
datesOut <- seq(1960, 2005, by = 5)

# NOTE: rule is an argument to approx(), which does the work.
# if not passed in, 1960 is returned as NAs (outside the range
# of dates given). In this case, rule = 2 assumes constant 
# equal to nearest neighbor, meaning 1960 data will be returned
# as equal to 1960.73 input data. So, function should only be 
# allowed to extrapolate for short distances.
interp(popmat, datesIn, datesOut, "linear", rule = 2)
interp(popmat, datesIn, datesOut, "exponential", rule = 2)
interp(popmat, datesIn, datesOut, "power", rule = 2)

# standard example data, taken from AGEINT PAS spreadsheet.
p1      <- c(100958, 466275, 624134, 559559, 446736, 370653, 301862, 249409,
  247473, 223014, 172260, 149338, 127242, 105715, 79614, 53660, 31021, 34596)
 p2    <- c(201916, 932550, 1248268, 1119118, 893472, 741306, 603724, 498818, 
 494946, 446028, 344520, 298676, 254484, 211430, 159228, 107320, 62042, 69192)

 # YYYY-MM-DD dates as character
 d1              <- "1980-04-07"
 d2              <- "1990-10-06"
 dout            <- "1985-07-01"
 
 (p_lin <- interp(cbind(p1,p2),c(d1,d2),dout,method = "linear"))
 (p_exp <- interp(cbind(p1,p2),c(d1,d2),dout,method = "exponential"))
 (p_pow <- interp(cbind(p1,p2),c(d1,d2),dout,method = "power"))
 
 # these were set to different dates: recreate above dates...
# expcheck <- c(198405.94950, 916338.81520, 1226568.46300, 1099663.56998,
# 877940.13607, 728419.34667, 593229.03315, 490146.68931, 486341.99907,
# 438274.37571, 338530.96200, 293483.90110, 250060.12230, 207754.56083,
# 156460.02559, 105454.37955, 60963.47946, 67989.18589)
#lincheck <- c(151279.0063,698682.8052,935224.2645,838462.8208,
#669404.8826,555399.4484,452320.6025,373723.1886,370822.2183,
#334171.9953,258120.4225,223773.2942,190663.8732,158407.0618,
#119296.4085,80406.0250,46482.9538,51839.8592)
lincheck - p_lin
# do controlled spreadsheet test, passing in decimal dates, etc

dd1 <- dec.date(d1)
# OK, this means we can do it all with Approx.
 
plot(EarlyPop, t ='l', col='red',
     ylim = c(400,1220000), 
     ylab = 'The counts', xlab = 'Age groups')

lines(interpolatePop(EarlyPop, LaterPop, EarlyDate, LaterDate, DesiredDate),
      t = 'l', col='black')

lines(LaterPop,
      t = 'l', col='blue')

legend(12,1000000, 
       legend = c('Pop in 1980-04-07',
                  'Pop in 1985-07-01',
                  'Pop in 1990-10-06'),
       col=c('red','black','blue'),
       lty = 1)
 
	

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
	if (grepl(method,pattern = "exp")){
		pattern <- "exponential"
	}
	if (grepl(method,pattern = "lin")){
		pattern <- "linear"
	}
	if (grepl(method,pattern = "pow")){
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


intlin <- interp(popmat = cbind(EarlyPop,LaterPop),
		datesIn = c(EarlyDate, LaterDate),
        datesOut = DesiredDate,
		method = "lin")
intexp <- interp(popmat = cbind(EarlyPop,LaterPop),
		datesIn = c(EarlyDate, LaterDate),
		datesOut = DesiredDate,
		method = "Exponential")
intpow2 <- interp(popmat = cbind(EarlyPop,LaterPop),
		datesIn = c(EarlyDate, LaterDate),
		datesOut = DesiredDate,
		method = "power")
intpow1.5 <- interp(popmat = cbind(EarlyPop,LaterPop),
		datesIn = c(EarlyDate, LaterDate),
		datesOut = DesiredDate,
		method = "power",power = 1.5)

intlin - intexp
plot(intlin,type='l')
lines(intpow2,col = "red")
lines(intexp,col = "blue")
lines(intpow1.5,col = "magenta")
legend("topright", lty = 1, col = c("black","red","blue"),
		legend = c("linear","power 2","exponential"))

# and again with matrices
LIN <- interp(popmat = cbind(p1,p2,p3),
		datesIn = c("1980-04-07", "1990-10-06",  "1997-10-01"),
		datesOut <- paste0(1980:1997,"-07-01"),
		method = "Linear")
EXP <- interp(popmat = cbind(p1,p2,p3),
		datesIn = c("1980-04-07", "1990-10-06",  "1997-10-01"),
		datesOut <- paste0(1980:1997,"-07-01"),
		method = "Exponential")
POW <- interp(popmat = cbind(p1,p2,p3),
		datesIn = c("1980-04-07", "1990-10-06",  "1997-10-01"),
		datesOut <- paste0(1980:1997,"-07-01"),
		method = "Power")

matplot(LIN - EXP, type = 'l')
matplot(LIN - POW, type = 'l')
matplot(EXP - POW, type = 'l')

for (i in 1:ncol(LIN)){
	plot(LIN[,i],type = 'l')
	lines(POW[,i],col = "red")
	lines(EXP[,i],col = "blue")
	locator(1)
}
