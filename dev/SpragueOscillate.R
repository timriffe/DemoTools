
# Author: tim
###############################################################################


# want some data w age heaping here:
Value <- c(9544406,7471790,11590109,11881844,11872503,12968350,11993151,10033918,
14312222,8111523,15311047,6861510,13305117,7454575,9015381,10325432,
9055588,5519173,12546779,4784102,13365429,4630254,9595545,4727963,
5195032,15061479,5467392,4011539,8033850,1972327,17396266,1647397,
6539557,2233521,2101024,16768198,3211834,1923169,4472854,
1182245,15874081,1017752,3673865,1247304,1029243,12619050,1499847,
1250321,2862148,723195,12396632,733501,2186678,777379,810700,
7298270,1116032,650402,1465209,411834,9478824,429296,1190060,
446290,362767,4998209,388753,334629,593906,178133,
4560342,179460,481230,159087,155831,1606147,166763,93569,182238,
53567,1715697,127486,150782,52332,48664,456387,46978,34448,
44015,19172,329149,48004,28574,9200,7003,75195,13140,5889,
18915,21221,72373)
Age <- 0:100
names(Value) <- Age
#barplot(Value, main = "yup, these have heaping!")
spragueOscillate(Value, Age, TRUE)
spragueSimple( groupAges(Value, Age))

spragueOscillate <- function(Value, Age, OAG = TRUE){
	
	N     <- length(Value)
	if (OAG){
		open   <- Value[N]
		OA     <- Age[N]
		Value  <- Value[-N]
		Age    <- Age[-N]
		N      <- N - 1
	} 
	
# select which ages to keep:
	p1x1   <- matrix(nrow = length(Value), ncol = 5)
	rownames(p1x1) <- Age
	for (i in 0:4){
		Age.i.5             <- calcAgeN(Age, shiftdown = i)
		keep.i              <- rep(rle(Age.i.5)$leng, rle(Age.i.5)$leng) == 5
		Age.i.5             <- Age.i.5[keep.i]
		Val.i               <- Value[keep.i]
		Val.i.5             <- groupAges(Val.i, Age.i, AgeN = Age.i.5)
		# make fake open age
		Val.i.5             <- c(Val.i.5, pi)
		names(Val.i.5)      <- c(unique(Age.i.5), max(Age.i.5) + 5)
		pop.est             <- spragueSimple(Val.i.5)
		pop.est[pop.est < 0] <- NA
		pop.est             <- pop.est[-length(pop.est)]
		p1x1[keep.i, i + 1] <- pop.est
	}
	
	p.out <- rowMeans(p1x1, na.rm = TRUE)
	if (OAG){
		Age   <- c(Age, OA)
		p.out <- c(p.out, open)
		names(p.out) <- Age
	}
	
	p.out
}

#spragueOscillate(Value, Age, TRUE)
#spragueSimple( groupAges(Value, Age))
#p1x12 <- p1x1
#p1x12[p1x12 < 0] <- NA
#names(Value) <- Age
#Val5 <- groupAges(Value, Age)
#matplot(p1x1,type='l', lty =1, col = gray(.5))
#lines(Value,  col = gray(.8))
#lines(grabill(as.matrix(Val5)),col = "green")
#lines(rowMeans(p1x1,na.rm=TRUE),col = "red",lwd=2)
##
#names(Value) <- Age
#Val5        <- groupAges(Value, Age)
#grabill(as.matrix(Val5))
#rowMeans(p1x12,na.rm=TRUE)
#write.csv(p1x1,file = "/home/tim/git/DemoTools/devdata/p1x1.csv")


