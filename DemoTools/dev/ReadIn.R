
# Author: tim
###############################################################################

setwd("/home/tim/workspace/DemoTools")

#D1 <- read.csv("data/Deaths_1x1.csv")
D5 <- read.csv("data/Deaths_5x1.csv")

dim(D5)
unique(D1$AgeSpan)
unique(D5$AgeSpan)
head(D5[D5$TimeLabel == 1817 & D5$IndicatorID == 195 & D5$SexID == 1, ])

D1 <- D5[D5$TimeLabel == 1820 & D5$IndicatorID == 195 & D5$SexID == 1 & D5$SubGroupTypeID == 2, ]
D1 <- D1[order(D1$AgeStart), ]

# D1$DataValue is the value column, deaths or population...
# AgeStart is the lower bound of the Age group
# AgeSpan is the interval of the age group (-1 is open)

getYear <- function(X){
# TimeLabel ought to be integer, and this should be the year,
# but not necessarily so. It is contained in the TimeStart column,
# which ought to be coercible to Date class, but which is coded
# with an ambiguous date. The strategy is therefore to take
# TimeMid and subtract half of TimeDuration. 
# This will not yield an integer if data aren't annualized
	
	X$Year <- X$TimeMid - X$TimeDuration / 2
	X
}
#D1 <- getYear(D1)
#head(D1)


# get an example

