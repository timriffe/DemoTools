# Author: tim
###############################################################################

# -----------------------------------------------------
#
blogit <- function(lx){
	lx <- lx/lx[1]
	Yx <- log((1-lx)/lx)/2
	Yx
}
bexpit <- function(Yx){
	1 / (1+exp(2*Yx))
}

HMD <- local(get(load("/home/tim/git/DistributionTTD/DistributionTTD/Data/HMDresults.Rdata")))
head(HMD)

lxs <- HMD$lx[HMD$CNTRY == "SWE" & HMD$Year == 2000 & HMD$Sex == "f"]/1e5

lxi <- HMD$lx[HMD$CNTRY == "JPN" & HMD$Year == 2010 & HMD$Sex == "f"]/1e5

Yxs <- blogit(lxs)[-1]
Yxi <- blogit(lxi)[-1]

plot(Yxs,Yxi)

abline(lm(Yxi~Yxs))
mod <- lm(Yxi~Yxs)

lxa <- c(1,bexpit(predict(mod,data.frame(Yxs=Yxs))))

a <- 0:110
plot(a,lxs,type='l')
lines(a, lxi,col="red")
lines(a, lxa, col = "blue",lty=2)

brass_adjust <- function(lx_in, age_in, lx_standard, age_standard=0:110){
	
	stopifnot(all(age_in %in% age_standard))
	
	# make ages match- no need for AgeInt if working only with lx...
	ind     <- age_standard %in% age_in
	lxs     <- lx_standard[ind] 
	
	# assuming first values 1, turning to Inf
	Yxs     <- blogit(lxs)[-1]
	Yxi     <- blogit(lx_in)[-1]
	
	mod     <- lm(Yxi~Yxs)
	
	# predict again for each age of the standard
	Yxi_out <- predict(mod,data.frame(Yxs=Yxs))
	lx_out  <- c(1,bexpit(Yxi_out))
	
	lx_out
}

plot(Yxs,Yxi)
abline(v=0);abline(h=0)
abline(a=0,b=1)
abline(lm(Yxi~Yxs))

