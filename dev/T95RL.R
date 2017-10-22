
# Author: Juan Galeano & Tim Riffe
###############################################################################

# this is here to preserve the old version / compare
T9R5L_old <- function(Value, Age, ns = 0){
	
	# Get Px and Pxsum (Px+)  
	PxPxsum<- function(x){
		
		Px    <- x[seq(1, length(x), 5)]
		A2    <- Px[ -length(Px)]
		sum04 <- as.numeric(tapply( x, (seq_along(x)-1) %/% 5, sum))
		sum04 <- sum04[-length(sum04)]
		Pxsum <- sum04 - A2
		aj1   <- list(A2, Pxsum, Px)
		aj1
	}
	
	inicio<-PxPxsum(Value) # create a list with the three elements 
	sumPxPxsum<-inicio[[1]]+inicio[[2]]
	Px<-inicio[[3]]
	
	#B<-inicio[[2]]
	#vector2<-c(1:(length(B)-1))  
	
	#for (i in c(1:(length(B)-1))) {
	#  vector2[i] <- (B[i]+B[i+1])
	#  vector2
	#}
	#vector2
	
	#calculate correction factor (CF) and P'x (Pxp) and P'x+ (Pxsump)
	
	f_AJUSTE <- function(A,B){
		
		CF<-8/9*((A+c(1,zoo::rollapply(B, 2, sum)))/c(1,zoo::rollapply(B, 2, sum)))
		CF[1]<-1
		Pxp<-c(A-(CF-1)*c(1,zoo::rollapply(B, 2, sum)))
		Pxsump<-(c(zoo::rollapply(CF, 2, sum),(CF[length(CF)]+1))-1)*B
		
		FACTORC <- CF
		POB1 <- Pxp
		POB2 <- Pxsump
		
		aj <- list(FACTORC,POB1,POB2)
		aj
	}
	
	# Iterate 
	Pxp    <- inicio[[1]]
	
	Pxsump <- inicio[[2]]
	

	# TR?????????????????????
	for (i in 1:100) {
		
		ajuste <- f_AJUSTE(Pxp,Pxsump)
		Pxp    <- ajuste[[2]] 
		Pxsump <- ajuste[[3]] 
		
		#i <- i + 1
	}
	
	finalc <- ajuste[[1]] 
	finalA <- ajuste[[2]]
	
	result <- list(finalc,finalA,Pxsump)
	
	df1 <-result[[2]]
	
	G<-(df1*.6)[c(2:(length(df1)-1))]
	H<-(df1*.4)[c(3:length(df1))]
	I<-G+H #f(x+2.5)
	
	Pxp5 <-c((inicio[[1]][1]+inicio[[2]][1]),I*5,utils::tail(sumPxPxsum,1),utils::tail(Px,1)) #5Pxp
	
	Px5<-c(Pxp5*(sum(sumPxPxsum)+utils::tail(Px,1)+ns)/sum(Pxp5))
	
	finalresult<-data.frame(Age_group=c(paste(seq(from = 0, to = length(Age), by = 5),
							seq(from = 4, to = length(Age)+4, by = 5), sep="-"),
					"Not Stated"),
			recorded=c(sumPxPxsum,utils::tail(Px,1),ns),
			corrected=c(Px5,NA))
	# TR: commented out this line because a not defined, broke code
	#finalresult$Age_group<-paste(as.character(a$Age_group))
	
	# TR: I commented this out because it was throwing NAMESPACE errors due to head and tail,
	# so maybe there's a more parismonious way to do this operation
	#finalresult[as.numeric(rownames(head(tail(finalresult[1],2),1))),1]<-paste(substr(finalresult[as.numeric(rownames(head(tail(finalresult[1],2),1))),1],1,2),"+",sep="")
	finalresult
}

