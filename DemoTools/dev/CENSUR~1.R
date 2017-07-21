
# Author: tim
###############################################################################

# based on Feeney's spreadsheet CENSUR~1.XLS
pop1 <- c(3831870,4502304,5397061,4630775,4193184,4114704,3770907,3274822,
2744786,2559755,2160716,1839025,1494043,1133409,870238,577972,313781,131547)
pop2 <- c(4292503,3988292,3852101,4492096,5347327,4571868,4190340,4085338,
3674127,3198934,2648360,2382691,1970485,1584699,1172155,736258,408191,53116)
Age  <- seq(0,85,by=5)
date1 <- "1960-10-01"
date2 <- "1970-10-01"

survRatio <- function(pop1, pop2, Age, date1, date2){
	stopifnot(all.equal(length(pop1),length(pop2),length(Age)))
	
	# census spacing used to stagger ages. What if censuses not exactly
	# 5 or 10 years apart?
	date1 <- dec.date(date1)
	date2 <- dec.date(date2)
	
}