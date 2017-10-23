# Author: Juan Galeano
###############################################################################

#' Feeney technique for removing "zigzag" from age data.

#' @description  Fenney technique for removing "zigzag" from age data. This comes from 
#' Feeney, G. 2013 "Removing "Zigzag" from Age Data," http://demographer.com/white-papers/2013-removing-zigzag-from-age-data/


#' @param Values numeric. A vector of demographic counts in five years age groups.
#' @param age character. A vector with ages (five years grups).
#' @param start numeric. The position in the age vector of the age group at the starting point of the zigzag.
#' @param end numeric. The position in the age vector of the age group at the ending point of the zigzag.

#' @details Five year age groups are assumed.

#' @return a dataframe with with the smoothed population distribution and the Weighted Sum of Squared Deviations . 
#' @export

#' @examples 
#' #data from feeney2013, Feeney, G. 2013 "Removing "Zigzag" from Age Data," http://demographer.com/white-papers/2013-removing-zigzag-from-age-data/

# age <-c("0","1-4","5-9","10-14","15-19","20-24","25-29","30-34", "35-39","40-44","45-49",
#        "50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+")
# Value <- c(13331, 4151,1746,1585,3859,8354,11146,12076,12216,12016,12473,
#            11513,12899,11413,12710,11516, 11408,6733,4031,2069)
#	
Value <- c(13331,4151,1746,1585,3859,8354,11146,12076,
	12216,12016,12473,11513,12899,11413,12710,11516,11408,6733,4031,2069)
Age <- c(0,1,seq(5,90,by=5))
	
 zigzag(Value, age,10,18)

 # function for average adjacent values excluding index value
avg_adj <- function(x){
	n <- length(x)
	(shift.vector(x, -1, NA) + shift.vector(x, 1, NA)) / 2
}

#avg_adj(Value)
#p <- c(0.0235802695087692,0.0724286618207911,0.0242327829742267,0.0883411499065237)
zigzag_min <- function(Value, Age, start = 40, end = 80, p = rep(0.05, (end - start) / 10)){
  
   # first, we need an odd number of age groups to smooth over.
   # we enforce this
	if ((((end - start)+5)/5) %% 2 == 0){
		if (start > 35){
			start <- start - 5
		} else {
			start <- start + 5
		}
	}
	
  # 1) this can be a multiplier later
  id <- rep(0,length(Age))
  id[Age >= start & Age <= end] <- 1
  
  Smoothed <- zigzagp(Value = Value, Age = Age, start = start, end = end, p = p)
  
  # 7) again adjacent avg
  avg2  <- avg_adj(Smoothed) * id
  avg2[is.na(avg2)] <- 0
  
  # 8) difference squared
  diffsq <- ((Smoothed * id) - avg2)^2

  # 9) return weighted sum sq dev
  sum(diffsq * 1 / Smoothed)
  
}

zigzagp <- function(Value, Age, start = 40, end = 80, p = rep(0.05, (end - start) / 10)){
	# first, we need an odd number of age groups to smooth over.
	# we enforce this
	if ((((end - start)+5)/5) %% 2 == 0){
		if (start > 35){
			start <- start - 5
		} else {
			start <- start + 5
		}
	}
	
	# 1) this can be a multiplier later
	id          <- rep(0,length(Age))
	id[Age >= start & Age <= end] <- 1
	
	# 2) average of value above and below
	avg1        <- avg_adj(Value) * id
	
	# 3) take difference
	diff1       <- (Value - avg) * id
	
	# 4) p used differently
	pid         <- Age %in% as.integer(seq((start + 5),end-5,by=10))
	npid        <- as.logical(id * !pid)
	npfromabove <- npid & Age < end
	npfrombelow <- npid & Age > start
	
	# 5) get Np
	Np          <- Value[pid] * p # show be 4 numbers in example
	
	# 6) use Np to smooth
	Smoothed                <- Value
	Smoothed[pid]           <- Value[pid] - Np
	Smoothed[npfromabove]   <- Smoothed[npfromabove] + Np / 2
	Smoothed[npfrombelow]   <- Smoothed[npfrombelow] + Np / 2

	# return smoothed
	Smoothed
}

zigzag <- function(Value, Age, start = 40, end = 80){
	p <- rep(0.05, (end - start) / 10)
	p <- optim(p , zigzag_min, Value = Value, Age = Age, start = start, end = end)$par
	Smoothed <- zigzagp(Value = Value, Age = Age, start = start, end = end, p = p)
	Smoothed
}


  #
  # STEP 1: CREATE A DF WITH AGES AND POPULATION COUNTS
  # DF<-data.frame(Age, Value)
  # 
  # # STEP 2: CREATE AND INDEX(P) OF THE SAME LENGTH OF THE AGE VECTOR
  # DF$P<-seq(1:length(DF$age))
  # 
  # # STEP 3: COMPUTE THE AVERAGE ADJACENT (AA) VALUES
  # DF$AA<-0
  # for(i in 2:length(DF$Value)){
  #   DF$AA[i]<-(DF$Values[i-1]+DF$Value[i+1])/2
  # }
  # 
  # # STEP 4: COMPUTE THE DIFFERENCE (diff) between the population values and the AA values
  # DF$diff<-with(DF,Value-AA)
  # 
  # # STEP 5: CREATE THE i COLUMN SETTING THE START AND END POSITION OF AGES WHERE WE NEED TO REMOVE THE ZIGZAG
  # DF$i<-with(DF, c(rep(NA, start-1),seq(1:(end-start+1)),rep(NA,length(age)-end)))
  # 
  # # STEP 6: CREATE A VECTOR WITH A FIXED VALUE FOR PARS (p) ON EVEN NUMBERS
  # DF$p<-with(DF, ifelse(i%%2==0,p,NA))
  # 
  # # STEP 7: COMPUTE Np
  # DF$Np<-with(DF,Values*p)
  # 
  # # STEP 8: CALCULATE THE SMOOTHED DISTRIBUTION 
  # DF$SMOOTHED<-with(DF, ifelse(is.na(i),Value,
  #                       ifelse(i==1,Value+sum(as.numeric(na.omit(DF[DF$i==2,"Np"]))/2),
  #                       ifelse(i%%2==0,Value-Np,
  #                       ifelse(i==3,  Value+sum(as.numeric(na.omit(DF[DF$i==2|DF$i==4,"Np"]))/2),
  #                       ifelse(i==5,  Value+sum(as.numeric(na.omit(DF[DF$i==4|DF$i==6,"Np"]))/2),
  #                       ifelse(i==7,  Value+sum(as.numeric(na.omit(DF[DF$i==6|DF$i==8,"Np"]))/2),
  #                       ifelse(i==9,  Value+sum(as.numeric(na.omit(DF[DF$i==8|DF$i==10,"Np"]))/2),
  #                       ifelse(i==11, Value+sum(as.numeric(na.omit(DF[DF$i==10|DF$i==12,"Np"]))/2),
  #                       ifelse(i==13, Value+sum(as.numeric(na.omit(DF[DF$i==12|DF$i==14,"Np"]))/2)
  #                                                                                      ,0))))))))))
  # 
  # # STEP 9: COMPUTE THE AVERAGE ADJACENT (AA2) VALUES FOR THE SMOOTHED DISTRIBUTION
  # AA2<-0
  # for(i in 2:length(DF$SMOOTHED)){
  #   
  #   AA2[i]<-(DF$SMOOTHED[i-1]+DF$SMOOTHED[i+1])/2
  #   
  # }
  # 
  # DF$AA2<-with(DF, ifelse(is.na(i), NA,AA2))
  # 
  # # STEP 10: COMPUTE THE SQUARED DIFFERENCE BETWEEN THE VALUES OF THE SMOOTHED DISTRIBUTION AND AA2
  # DF$SDIFF<-with(DF, (SMOOTHED-AA2)^2)
  # 
  # # STEP 11: COMPUTE THE Weighted Sum of Squared Deviations (WSSD)
  # WSSD<-sum(1/DF$SMOOTHED[start:end]*na.omit(DF$SDIFF))
  # 
  # #STEP 12: RETURN THE COMPLETE DATAFRAME AND THE VALUE OF THE WSSD
  # print(DF)
  # print(paste("Weighted Sum of Squared Deviations", round(WSSD,0),sep=": "))

