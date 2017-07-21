
# Author: tim
###############################################################################

# based on Feeney's spreadsheet CENSUR~1.XLS

survRatio <- function(pop1, pop2, Age1, Age2, date1, date2){
	
	# census spacing used to stagger ages. What if censuses not exactly
	# 5 or 10 years apart?
	date1 <- dec.date(date1)
	date2 <- dec.date(date2)
	
}