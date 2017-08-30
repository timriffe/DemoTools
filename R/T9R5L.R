# Author: Juan Galeano
###############################################################################

#' Feeney T9R5L formula on 9 years to correct for heaping on multiples of five

#' @description  Fenney technique for correcting age distributions for heaping on multiples of five. This comes from 
#' Feeney, G. 1979. "A technique for correcting age distributions for heaping on multiples of five," 
#' Paper presented at Asian and Pacific Census Forum. Vol. 5:12-14
 
#' @param Value numeric. A vector of demographic counts in single age groups.
#' @param Age numeric or character. A vector with ages (in single years).
#' @param ns numeric. In case their is some Age not stated population. By default this is equal to 0.

#' @details Single year age groups are assumed.

#' @return a dataframe with 3 columns: Five year age groups, recoded values and corrected values. 
#' @export

#' @references 
#' \insertRef{feeney1979}{DemoTools}

#' @examples 
#' # data from gray1987missingages, Table 2, page 21: Aplication of Q1 and Q2 linear operators, Bangladesh Census, 1 March 1974-Males.
#' Pop <-c(2337,3873,3882,3952,4056,3685,3687,3683,3611,3175,
#'         3457,2379,3023,2375,2316,2586,2014,2123,2584,1475,
#'         3006,1299,1236,1052,992,3550,1334,1314,1337,942,
#'         3951,1128,1108,727,610,3919,1221,868,979,637,
#'         3409,887,687,533,313,2488,677,426,524,333,
#'         2259,551,363,290,226,1153,379,217,223,152,
#'         1500,319,175,143,89,670,149,96,97,69,
#'         696,170,60,38,23,745,15)
#' Ages<-c(0:74,"75+","Not Stated")
#' T9R5L(Pop, Ages,15)
T9R5L<- function(Value,Age,ns=0){
  
  # Get Px and Pxsum (Px+)  
  PxPxsum<- function(x){
    library(zoo)
    Px <- x[seq(1, length(x), 5)]
    A2<-Px[-length(Px)]
    sum04<-as.numeric(tapply( x, (seq_along(x)-1) %/% 5, sum))
    sum04<-sum04[-length(sum04)]
    Pxsum<-sum04-A2
    aj1 <- list(A2,Pxsum, Px)
    aj1}
  
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
    
    CF<-8/9*((A+c(1,rollapply(B, 2, sum)))/c(1,rollapply(B, 2, sum)))
    CF[1]<-1
    Pxp<-c(A-(CF-1)*c(1,rollapply(B, 2, sum)))
    Pxsump<-(c(rollapply(CF, 2, sum),(CF[length(CF)]+1))-1)*B
    
    FACTORC <- CF
    POB1 <- Pxp
    POB2 <- Pxsump
    
    aj <- list(FACTORC,POB1,POB2)
    aj
  }
  
  # Iterate 
  Pxp<-inicio[[1]]
  
  Pxsump<-inicio[[2]]
  
  i <- 1
  
  for (i in 1:100) {
    
    ajuste <- f_AJUSTE(Pxp,Pxsump)
    Pxp <- ajuste[[2]] 
    Pxsump <- ajuste[[3]] 
    
    i <- i + 1
  }
  
  finalc <- ajuste[[1]] 
  finalA<- ajuste[[2]]
  
  result <- list(finalc,finalA,Pxsump)
  
  df1 <-result[[2]]
  
  G<-(df1*.6)[c(2:(length(df1)-1))]
  H<-(df1*.4)[c(3:length(df1))]
  I<-G+H #f(x+2.5)
  
  Pxp5 <-c((inicio[[1]][1]+inicio[[2]][1]),I*5,tail(sumPxPxsum,1),tail(Px,1)) #5Pxp
  Px5<-c(Pxp5*(sum(sumPxPxsum)+tail(Px,1)+ns)/sum(Pxp5))
  
  finalresult<-data.frame(Age_group=c(paste(seq(from = 0, to = length(Age), by = 5),
                                            seq(from = 4, to = length(Age)+4, by = 5), sep="-"),"Not Stated"),
                          recorded=c(sumPxPxsum,tail(Px,1),ns),
                          corrected=c(Px5,NA))
  finalresult$Age_group<-paste(as.character(a$Age_group))
  finalresult[as.numeric(rownames(head(tail(finalresult[1],2),1))),1]<-paste(substr(finalresult[as.numeric(rownames(head(tail(finalresult[1],2),1))),1],1,2),"+",sep="")
  finalresult
}