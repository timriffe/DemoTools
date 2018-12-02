
# Author: tim
###############################################################################


#' Life expectancy in the open age group.
#' 
#' @description Get an estimate of life expectancy in the open age group.
#' @details This method estimates life expectancy in the open age group by fitting one of several potential old-age parametric mortality models, extrapolating rates to age 130, then backing out the implied remaining life expectancy in the open age group.
#' @inheritParams extra_mortality
#' @return life expectancy in the open age group
#' @references 
#' \insertRef{mortpak1988}{DemoTools}
#' @export
#' @examples
#' nMx <- c(0.12846,0.02477,0.00603,0.0034,
#'  		0.00417,0.00513,0.00581,0.00645,0.00725,
#'  		0.00813,0.00913,0.01199,0.01647,
#'  		0.0256,0.04047,0.06624,0.10638,0.19611)
#' Age <- c(0,1,seq(5,80,by=5))
#' 
#' 
#' aomegaMortalityLaws(nMx,Age,"Kannisto")
#' aomegaMortalityLaws(nMx,Age,"Kannisto_Makeham")
#' aomegaMortalityLaws(nMx,Age,"Makeham")
#' aomegaMortalityLaws(nMx,Age,"Gompertz")
#' aomegaMortalityLaws(nMx,Age,"GGompertz")
#' aomegaMortalityLaws(nMx,Age,extrapLaw="Beard")
#' aomegaMortalityLaws(nMx,Age,"Beard_Makeham")
#' aomegaMortalityLaws(nMx,Age,"Quadratic")

aomegaMortalityLaws <- function(
		mx, 
		Age, 
		law = c("kannisto",
				"kannisto_makeham",
				"gompertz", 
				"ggompertz", 
				"makeham", 
				"beard", 
				"beard_makeham", 
				"quadratic")[1],
		extrapFrom = max(Age),
		extrapFit = Age[Age >= 60],
		...){
	
	extrapLaw <- tolower(extrapLaw)
	OA        <- max(Age)
	x_extr    <- seq(OA, 130, by = .1)
	
	Mxnew <- extra_mortality(
			mx = nMx, 
			x = Age, 
			x_fit = extrapFit,
			x_extr = x_extr,
			law = extrapLaw,
			...)
	
	mx <- Mxnew$values[names2age(nMxext)>=OA]
	
	# acceptable approximation for small intervals
	lx <- c(1, exp(-cumsum(mx / 10)),0)
	# divide by 10 (interval) and 2 (avg), so 20.
	sum(shift.vector(lx, -1) + lx) / 20

}





