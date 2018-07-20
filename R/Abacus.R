# Author: tim
# a temp function family to make use of the Mortpak abacus abridged lifetable,
# mostly to get its closeout behavior.
###############################################################################

#' calculate simple linear model by hand
#' @description A legacy function internal to the mortpak lifetable. It calculates intercept, slope, correlation coef, and the residual sum of squares. Just like \code{lm()}? Unchecked.
#' @details Used internally by \code{AbacusLIFTB()}. 
#' @param X x
#' @param Y y
#' @param N length of X and Y
#' @return A list containing, A, B, R, SUMSQ, and potentially an error code.
#' @export
#' @references 
#' \insertRef{mortpak1988}{DemoTools}


REGRES <- function(X,Y,N)  { #  function(X,Y,N,A,B,R,SUMSQ,ErrorCode)
	#     -----------------------------------------------------------------------
	#     X,Y ARE INPUT DATA ARRAYS.
	#     N IS THE NUMBER OF INPUT DATA ELEMENTS (N <= length of X and Y). 
	
	#     A,B - FOR EQUATION Y=A*X+B
	#     R - COEFFICIENT OF CORRELATION.
	#     SUMSQ - SUM OF SQUARES.
	#     ErrorCode
	#       - 0 is no error
	#       - 1 is singular matrix (i.e. divide by 0)
	#       - 2 is N > lengtha(X)
	#       - 3 is N > length(Y)
	#       - 4 is N<=1, must be at least 2 data points 
	#     -----------------------------------------------------------------------
	
	# Initialize all output values, incase error condition
	A 		<- 0.0
	B 		<- 0.0
	R 		<- 0.0
	SUMSQ 	<- 0.0
	ErrorCode <- 0
	if(N>length(X)) {ErrorCode=2}
	if(N>length(Y)) {ErrorCode=3}
	if(N<=1) {ErrorCode=4}
	if(ErrorCode==0) {
		SUMX 	<- sum(X[1:N])
		SUMY 	<- sum(Y[1:N])
		SUMXX 	<- sum(X[1:N]^2)
		SUMYY 	<- sum(Y[1:N]^2)
		CROSP 	<- sum(X[1:N]*Y[1:N])
		SX 		<- (N*SUMXX-SUMX*SUMX)
		if(SX==0.0) {ErrorCode=1}
	}
	# If error condition observed, no further calculations. 
	if(ErrorCode==0) {
		SY 		<- (N*SUMYY-SUMY*SUMY)  
		A 		<- (N*CROSP-SUMX*SUMY)/SX
		B 		<- (SUMY*SUMXX-SUMX*CROSP)/SX
		R 		<- (N*CROSP-SUMX*SUMY)/sqrt(SX*SY)
		SUMSQ 	<- 0.0
		for (J in 1:N) {SUMSQ <- SUMSQ + (A*X[J]+B-Y[J])**2}
	}
	return(list("A"=A,"B"=B,"R"=R,"SUMSQ"=SUMSQ,"ErrorCode"=ErrorCode))
}

#' literal translation of mortpak from Fortran
#' @description This has not been modified in any meaningful way by the package maintainer, nor will it be documented in depth, as it is intended to be replaced in the future with a refactored version of the same. 
#' @details Try \code{AbacusLIFTB_wrap()} if you just want to make quick use of this functionality. Standard abridged ages are assumed. The lifetable stops at 100+. Note \code{QXMX} needs to be an object in single ages
#' @param NFIN integer the open age
#' @param NUMOUT the desired open age group (under 100)
#' @param NTYPE 1 means nqx. 2 means nMx. 3 means nMx but do not extrapolate to close out.
#' @param NSEX 1 is male. 2 is female
#' @param QXMX numeric vector of either nMx or nqx in abridged ages
#' @export 
#' 
#' @references 
#' \insertRef{mortpak1988}{DemoTools}

AbacusLIFTB <- function(NFIN,NUMOUT,NTYPE,NSEX,QXMX)  {
	# ---------------------------------------------------------------------------------------------------------------
	#         LIFTB is derived from MORTPAK software package and customized for Abacus
	#               UNITED NATIONS SOFTWARE PACKAGE FOR MORTALITY MEASUREMENT     
	# ---------------------------------------------------------------------------------------------------------------
	# NFIN - The open age group for the input q(x,n) or m(x,n) data.
	# NUMOUT - The open age group of the output table. Default to NFIN if value is 0.
	#
	#       Add extrapolation on mx if needed.
	# NTYPE - Type of input data where 1 is q(x,n) and 2 is m(x,n).  Values 3 is mx, but use open age group and do not extrapolate q.
	# NSEX -  Gender of input data where 1 is males and 2 is females.
	# QXMX -  Contains input m(x,n) or q(x,n) values identified by value of NTYPE. Vector dimensioned to 101. If NTYPE is 3, then
	#         QXMX(NFIN) must have a value at the open age group because it will not be estimaed by extrapolation q until no survivors remain.
	# 
	# ARRAY - Output life table dimensioned to ARRAY(101,9).  
	#         This format was used in original mortpak.  For return array, can use style suitable for projection program.
	#         ARRAY(age+1,colnum) where columns are m,q,l(x),d,L,s,T,E,a
	# 
	#     ErrorCode
	#       - 0 No error detected
	#       - 1 NFIN out of range of 65 to 100 (must be multiple of 5)
	#       - 2 NTYPE is not 1 or 2
	#       - 3 NSEX is not 1 or 2
	#       - 4 Input q(x) or m(x) is outside the range of 0 and 1
	#       - 5 Warning - convergence not reached when computing regression in completion of life table 
	#       - 6 Warning - mortality rates not increasing monotonically at older ages
	# -----------------------------------------------------------------------------------------------------------------
	
	ARRAY 		<- array(rep(0.0,909),dim=c(101,9))  # Dimension and initialize output array.
	QVAL			<- rep(0.0,25)
	XVAL			<- rep(0.0,6)
	TT			<- rep(0.0,3)
	AM			<- rep(0.0,101)
	A				<- rep(0.0,101)
	G				<- rep(0.0,101)
	XV			<- rep(0.0,6)
	YV			<- rep(0.0,6)
	oag			<- rep(0.0,9)
	XAgeGroup 	<- rep(0.0,22)
	YLnMxValues 	<- rep(0.0,22)
	
	# New feature, when NTYPE is 3, do not extrapolate qx but instead use mx at open age group (e.g. ax(oag)<-1.0/mx(oag))
	UseMxOpenAgeGroup <- 0   # Extrapolate q until no survivors remain.  The mx at the open age group is based on extrapolated qx.  
	if(NTYPE == 3)  UseMxOpenAgeGroup <- 1  #  Instead of extrapolating, use open age group from mx
	if(NTYPE == 3)  NTYPE <- 2    # Set NTYPE to traditional value for mx
	
	
	# Initialize values - when iteration starts, program stops when not defined for example   C3 <- C2
	C1 	<- 0.0
	C2 	<- 0.0
	C3 	<- 0.0
	SUM1 	<- 0.0
	SUM2 	<- 0.0
	SUM3 	<- 0.0
	
	
	nfnout 	<- NUMOUT
	if(NUMOUT == 0) nfnout <- NFIN
	NEND 		<- NFIN+1
	ErrorCode	<- 0
	QXMX[3]	<- 0.0
	QXMX[4]	<- 0.0
	QXMX[5]	<- 0.0
	ARRAY[100,2]<-1.0
	
	if(NFIN < 65 | NFIN > 100) ErrorCode <- 1  #  Last age group must be between 60-65 and 95-100.
	#  If NTYPE is originally set to 3 it will have a value ot 2 here.  However variable UseMxOpenAgeGroup will be set to 1.
	if(NTYPE < 1 | NTYPE > 2)  ErrorCode <- 2   #  TYPE OF MORTALITY RATE MUST BE A 1 TO 2.  
	if(NSEX < 1  | NSEX > 2)   ErrorCode <- 3    #  sex code MUST BE A 1 TO 2.
	
	ErrFlag<-0  #  Note that error codes 1-4 are errors which stops the function while 5 and 6 are warnings (i.e. execution continues)
	for (I in 1:2) {
		if(QXMX[I] <= 0.0 | QXMX[I] >= 1.0) ErrFlag<-1
	}
	#QXMX[seq(6,NFIN, by=5)]
	# Qx needs to be between 0 and 1
	for (I in seq(6,NFIN, by=5)) {
		if(QXMX[I] <= 0.0 | QXMX[I] >= 1.0) ErrFlag<-1
	}
	if(ErrFlag != 0) {
		ErrorCode <- 4
	}  #     ERROR IN LIFTB.  INPUT Q(X,N) OR M(X,N) VALUE(S) IS OUTSIDE THE RANGE OF ZERO AND ONE.  
	if (ErrorCode != 0) {
		return(list("ARRAY"=ARRAY,"ErrorCode"=ErrorCode))
	}   # 
	
	if (NTYPE == 1) {
		for (I in seq(6,NEND,by=5))  {ARRAY[I,2]<-QXMX[I]}
		ARRAY[1:5,2]	<- QXMX[1:5]
		ARRAY[NEND,2]	<- 1.0
		A[6]			<- 2.5
		A[11]			<- 2.5
		NM10			<- NEND-10
		NM5				<- NEND-5
		for (I in seq(16,NM5,by=5)) {
			A[I]<-2.5
		}
		for (IT in 1:20)     {
			for(I in seq(11,NM5,by=5))  {
				AM[I]=QXMX[I]/(5.0-(5.0-A[I])*QXMX[I])
			}
			for (I in seq(16,NM10,by=5))  {
				AK		<-log(AM[I+5]/AM[I-5])/10.0
				A[I]	<-2.5-(25.0/12.0)*(AM[I]-AK)
				if(A[I] < 1.0) A[I] <- 1.0
				if(I >= 61) {
					tmp 					<- 0.8*A[I-5]
					if(A[I] < tmp) A[I]	<- tmp
				}
			}
			A[NM5] 					<- 2.5-(25.0/12.0)*(AM[NM5]-AK)
			if(A[NM5] < 1.0) A[NM5] 	<- 1.0
			tmp 						<- 0.8*A[NM5-5]
			if(A[NM5] < tmp) A[NM5] 	<- tmp
		}
	}
	if(NTYPE == 2)  {
		A[6]		<- 2.5
		A[11]		<- 2.5
		NM10		<- NEND-10
		NM5			<- NEND-5
		for (I in seq(16,NM10,by=5)) {
			AK        <- log(QXMX[I+5]/QXMX[I-5])/10.0
			A[I]      <- 2.5-(25.0/12.0)*(QXMX[I]-AK)
		}
		#  Estimate ax, for the age group 5 years before the open age group.  To use the formula, we need to estimate 5m100.
		#  Estimate 5m100 by extrapolating ln(mx) and keep resultant value at or below m100+ (i.e. open age group)
		#  Problem with this method is that when 5M100 approaches one, it curves to the right.
		#  This might overstating 5M100 value and possibly be greater the input open age group M100+
		if(UseMxOpenAgeGroup == 1)  {     
			#  Do regression on ln(mx)
			for (I in 1:6)  {  # For age groups 70 to 95
				AgeValue 		<- 65+5*I
				XAgeGroup[I] 	<- AgeValue
				YLnMxValues[I] 	<- log(QXMX[AgeValue+1])
			}
			tmp      <- REGRES(XAgeGroup,YLnMxValues,6)
			EstMx100 <- exp(tmp$A * 100 + tmp$B)   #  Regression on ln(mx), restore using exponent
			
			# Possibly bound value 5m95 and m100+
			#  Another possibility is to use a proportion of the open age group to estimate 5M100
			TmpVal   <- exp(8.0*QXMX[NM5+5]+0.2)
			EstMx100 <- QXMX[NM5+5]*TmpVal/(1.0+TmpVal)  # Logit function prevents values over 1.0
			
			# At this stage simply use mx from the open age group since at age group 100, the 5 year
			# and open age group are close together.  Afterwards, use a better method.  For example
			#  multiply the mx at the open age group by .92 to .98 as a better estimate.
			
			#   AK<-log(EstMx100/QXMX[NM5-5])/10.0   
			AK<-log(QXMX[NM5+5]/QXMX[NM5-5])/10.0   # This uses open age group, in many cases seems to match slon results.
			
			
		}
		
		A[NM5]<-2.5-(25.0/12.0)*(QXMX[NM5]-AK)   # For MORTPAK, use AK from previous age group to avoid open age group
		
		for (I in seq(6,NFIN,by=5)) {
			ARRAY[I,2]<-5.0*QXMX[I]/(1.0+(5.0-A[I])*QXMX[I])
		}
		ARRAY[NEND,2]<-1.0
		if(NSEX == 1 & QXMX[1] < 0.1072) {
			ARRAY[1,2]<-((1.0+QXMX[1]-.0425*QXMX[1])-sqrt((.0425*QXMX[1]-QXMX[1]-1.0)*(.0425*QXMX[1]-QXMX[1]-1.0)-4.0*2.875*QXMX[1]*QXMX[1]))/(2.0*2.875*QXMX[1])
			ARRAY[2,2]<-4.0*QXMX[2]/(1.0+(4.0-1.653+3.013*ARRAY[1,2])*QXMX[2])
		}
		if(NSEX == 1 & QXMX[1] >= 0.1072)  {
			ARRAY[1,2]<-1.0*QXMX[1]/(1.0+.67*QXMX[1])
			ARRAY[2,2]<-4.0*QXMX[2]/(1.0+2.648*QXMX[2])
		}
		if(NSEX == 2 & QXMX[1] < 0.1072) {
			ARRAY[1,2]<-((1.0+QXMX[1]-.0500*QXMX[1])-sqrt((.0500*QXMX[1]-QXMX[1]-1.0)*(.0500*QXMX[1]-QXMX[1]-1.0)-4.0*3.000*QXMX[1]*QXMX[1]))/(2.0*3.000*QXMX[1])
			ARRAY[2,2]<-4.0*QXMX[2]/(1.0+(4.0-1.524+1.627*ARRAY[1,2])*QXMX[2])
		}
		if(NSEX == 2 & QXMX[1] >= 0.1072) {
			ARRAY[1,2]<-1.0*QXMX[1]/(1.0+.65*QXMX[1])
			ARRAY[2,2]<-4.0*QXMX[2]/(1.0+2.639*QXMX[2])
		}
		ARRAY[3,2]=0.0
		ARRAY[4,2]=0.0
		ARRAY[5,2]=0.0
	}
	ARRAY[1,3]<-100000.0
	for (I in 2:6) {
		ARRAY[I,3]<-ARRAY[I-1,3]*(1.0-ARRAY[I-1,2])
		ARRAY[I-1,4]=ARRAY[I-1,3]*ARRAY[I-1,2]
	}
	for (I in seq(11,NEND,by=5)) {
		ARRAY[I,3]<-ARRAY[I-5,3]*(1.0-ARRAY[I-5,2])
		ARRAY[I-5,4]<-ARRAY[I-5,3]*ARRAY[I-5,2]
	}
	ARRAY[NEND,4]<-ARRAY[NEND,3]
	if(ARRAY[1,2] >= 0.100 & NSEX == 1) {
		ARRAY[1,5]<-0.33*ARRAY[1,3]+0.67*ARRAY[2,3]
		ARRAY[2,5]<-1.352*ARRAY[2,3]+2.648*ARRAY[6,3]
	}
	if(ARRAY[1,2] < 0.100 & NSEX == 1) {
		ARRAY[1,5]<-(0.0425+2.875*ARRAY[1,2])*ARRAY[1,3]+(0.9575-2.875*ARRAY[1,2])*ARRAY[2,3]
		ARRAY[2,5]<-(1.653-3.013*ARRAY[1,2])*ARRAY[2,3]+(2.347+3.013*ARRAY[1,2])*ARRAY[6,3]
	}
	if(ARRAY[1,2] >= 0.100 & NSEX == 2) {
		ARRAY[1,5]<-0.35*ARRAY[1,3]+0.65*ARRAY[2,3]
		ARRAY[2,5]<-1.361*ARRAY[2,3]+2.639*ARRAY[6,3]
	}
	if(ARRAY[1,2] < 0.100 & NSEX == 2) {
		ARRAY[1,5]<-(0.0500+3.000*ARRAY[1,2])*ARRAY[1,3]+(0.9500-3.000*ARRAY[1,2])*ARRAY[2,3]
		ARRAY[2,5]<-(1.524-1.627*ARRAY[1,2])*ARRAY[2,3]+(2.476+1.627*ARRAY[1,2])*ARRAY[6,3]
	}
	ARRAY[3,5]<-0.0
	ARRAY[4,5]<-0.0
	ARRAY[5,5]<-0.0
	for(I in seq(6,NFIN,by=5)) {ARRAY[I,5]<-A[I]*ARRAY[I,3]+(5.0-A[I])*ARRAY[I+5,3]}
	
	
	#  This block extrapolates qx to get Lx at the open age group. Skip if NTYPE is set to 3
	if(UseMxOpenAgeGroup == 0)  {     #  This block extrapolates qx to get Lx at the open age group. Skip if NTYPE is set to 3
		
		ErrorFlag<-0
		for (I in seq(46,NFIN,by=5)) {if(ARRAY[I,2] < ARRAY[I-5,2]) ErrorFlag<-1}
		if(ErrorFlag != 0)  ErrorCode<-6   # WARNING:  MORTALITY RATES AT OLDER AGES DO NOT INCREASE MONOTONICALLY. 
		NBEG<-NEND-30
		KN<-0
		for (I in seq(NBEG,NFIN,by=5)) {
			KN<-KN+1
			QVAL[KN]<-ARRAY[I,2]/(1.0-ARRAY[I,2])
			XVAL[KN]<-I-1
		}
		TMP1 <- (QVAL[5]-QVAL[3])/(QVAL[3]-QVAL[1])
		TMP2 <- (QVAL[6]-QVAL[4])/(QVAL[4]-QVAL[2])
		
		
		#  for (ii in 1:101) {
		#    cat("\n",spaces(3),sprintf("%2i  %10.5f%10.5f%10.0f%10.0f%10.0f%10.5f%10.0f%10.3f%10.3f",ii-1,ARRAY[ii,1],ARRAY[ii,2],ARRAY[ii,3],ARRAY[ii,4],ARRAY[ii,5],ARRAY[ii,6],ARRAY[ii,7],ARRAY[ii,8],ARRAY[ii,9]))
		
		#  }
		
		
		
		if(TMP1 <= 0  | TMP2 <= 0) {STVAL<-1.12} else {
			STVAL=(TMP1^0.1+TMP2^0.1)/2.0
		}
		#----------------------------------------------
		#     CALCULATION OF NON-LINEAR REGRESSION
		#----------------------------------------------
		DELTA<-0.05
		TT[3]<-STVAL-DELTA
		ErrorFlag<-1  #  If iteration converges, then set code to zero. 
		for (ITER in 1:99)  {
			for (I in 1:6)  {
				XV[I]<-TT[3]^XVAL[I]
				YV[I]<-QVAL[I]
			}
			FN<-6.0
			SUMX<-sum(XV[1:6])
			SUMY<-sum(YV[1:6])
			SUMXX<-sum(XV[1:6]^2)
			CROSP<-sum(XV[1:6]*YV[1:6])
			SX=(FN*SUMXX-SUMX*SUMX)
			if(SX == 0.0) break   #  If SX is zero, break to avoid divide by zero. Then set ErrorCode to 5 for not converging.
			TT[1]=(SUMY*SUMXX-SUMX*CROSP)/SX
			TT[2]=(FN*CROSP-SUMX*SUMY)/SX
			SUMSQ<-0.0
			for (J in 1:6)  {SUMSQ<-SUMSQ+(TT[1]+TT[2]*XV[J]-YV[J])^2}
			#  An epsilon value of TT(3)/100000.0 will produce enough accuracy. The epsilon value was reduced so that
			#  the regression will produce identical results to NLS, the original non-linear least squares procedure
			if(ITER >= 3)  {
				EPS<-TT[3]/500000.0  
				DIFF<-abs(TT[3]-C1)
				if(DIFF < EPS) {ErrorFlag<-0; break}
				DIFF<-abs(TT[3]-C2)
				if(DIFF < EPS) {ErrorFlag<-0; break}
			}
			C3<-C2
			C2<-C1
			C1<-TT[3]
			SUM3<-SUM2
			SUM2<-SUM1
			SUM1<-SUMSQ
			if(ITER == 1) TT[3]<-STVAL+DELTA
			if(ITER == 2) TT[3]<-STVAL
			if(ITER >= 3) {
				TT[3]<-((SUM2-SUM1)*(C3*C3-C2*C2)-(SUM3-SUM2)*(C2*C2-C1*C1))/2.0
				TT[3]<-TT[3]/((SUM2-SUM1)*(C3-C2)-(SUM3-SUM2)*(C2-C1))
				if(TT[3] < 0.50)TT[3]<-0.5
				if(TT[3] > 2.1) TT[3]<-2.1
			}
		}   #  Iteration loop
		if(TT[3] <= 1.0 | TT[3] >= 2.0)  ErrorFlag=1   #  Value for TT[3] out of range 1 to 2, give ErrorCode=5 for non-convergence
		if (ErrorFlag==1) {ARRAY[100,2]=0.0; ErrorCode<-5}   #  convergence not reached for completion of life table
		#---------------------------------------------------------------------
		#     COEFFICIENTS FOR NON-LINEAR REGRESSION EQUATION WERE CALCULATED.
		#---------------------------------------------------------------------
		KN<-KN+1
		KNP12<-KN+13
		IK<-NEND-6
		for (I in KN:KNP12)  {
			IK<-IK+5
			QVAL[I]<-TT[1]+TT[2]*TT[3]^IK
			QVAL[I]<-QVAL[I]/(1.0+QVAL[I])
		}
		if(NTYPE == 2) AM[KN-1]<-QXMX[NM5]
		if(NTYPE == 1) AM[KN-1]<-AM[NM5]
		AM[KNP12+1] <- 1.0
		G[KN:KNP12] <- 2.5
		for (IT in 1:10)  {
			for (I in KN:KNP12)  {
				if(G[I] <= 0.0) G[I]<-0.5
				AM[I]<-QVAL[I]/(5.0-(5.0-G[I])*QVAL[I])
			}
			for (I in KN:KNP12)  {
				AK<-log(AM[I+1]/AM[I-1])/10.0
				G[I]<-2.5-(25.0/12.0)*(AM[I]-AK)
			}
		}
		for (I in KN:KNP12)  {if(G[I] <= 0.0) G[I]<-0.5}
		CAPL<-0.0
		SMLX<-ARRAY[NEND,3]
		for (I in KN:KNP12) {
			SMLX1<-SMLX*(1.0-QVAL[I])
			#     CHECK SMLX1 FOR POSSIBLE UNDERFLOW IN EXTREME CASES.
			if(SMLX1 < 0.001) SMLX1<-0.0
			ARRAY[5*I-28,5]<-G[I]*SMLX+(5.0-G[I])*SMLX1
			ARRAY[5*I-28,1]<-AM[I]
			CAPL<-CAPL+G[I]*SMLX+(5.0-G[I])*SMLX1
			SMLX<-SMLX1
		}
		ARRAY[NEND,5]<-CAPL
		ARRAY[NEND,7]<-CAPL
		
		
		
	}  # End of block when UseMxOpenAgeGroup is 0.  This block extrapolates qx to get Lx at the open age group. Skip if NTYPE is set to 3.
	if(UseMxOpenAgeGroup == 1)  {
		CAPL <- ARRAY[NEND,3]/QXMX[NEND]
		ARRAY[NEND,5]<-CAPL
		ARRAY[NEND,7]<-CAPL
	}
	
	
	
	for (I in 1:5)  {
		if (ARRAY[I,5] >= 1.0) ARRAY[I,1] <-ARRAY[I,4]/ARRAY[I,5]
		if (ARRAY[I,5] < 1.0) ARRAY[I,1] <-0.0
	}
	for( I in seq(6,NEND,by=5)) {ARRAY[I,1]<-ARRAY[I,4]/ARRAY[I,5]}
	K<-NEND-5
	for (I in seq(6,K,by=5)) {
		J<-NEND-I+1
		ARRAY[J,7]=ARRAY[J+5,7]+ARRAY[J,5]
	}
	for(I in 1:5)  {
		J<-6-I
		ARRAY[J,7]<-ARRAY[J+1,7]+ARRAY[J,5]
	}
	SUM<-sum(ARRAY[1:5,5])
	ARRAY[1,6]<-SUM/500000.0
	ARRAY[2,6]<-ARRAY[6,5]/SUM
	KM10<-NEND-10
	for (I in seq(6,KM10,by=5))  {ARRAY[I,6]<-ARRAY[I+5,5]/ARRAY[I,5]}
	ARRAY[K,6]<-ARRAY[NEND,7]/ARRAY[K,7]
	ARRAY[1:5,8]<-ARRAY[1:5,7]/ARRAY[1:5,3]
	for(J in seq(6,NEND,by=5))  {ARRAY[J,8]<-ARRAY[J,7]/ARRAY[J,3]}
	for(J in seq(6,NFIN,by=5))  {ARRAY[J,9]<-A[J]}
	ARRAY[NEND,9]<-1.0/ARRAY[NEND,1]
	ARRAY[3,9]<-0.47
	ARRAY[4,9]<-0.49
	ARRAY[5,9]<-0.50
	if(ARRAY[1,2] >= 0.100 & NSEX == 1) {
		A[1] <- 0.33
		A[2] <- 1.352
	}
	if(ARRAY[1,2] >= 0.100 & NSEX == 2) {
		A[1] <- 0.35
		A[2] <- 1.361
	}
	if(ARRAY[1,2] < 0.100 & NSEX == 1) {
		A[1] <- 0.0425+2.875*ARRAY[1,2]
		A[2] <- 1.653-3.013*ARRAY[1,2]
	}
	if(ARRAY[1,2] < 0.100 & NSEX == 2) {
		A[1] <- 0.050+3.00*ARRAY[1,2]
		A[2] <- 1.524-1.627*ARRAY[1,2]
	}
	ARRAY[1,9]<-A[1]
	ARRAY[2,9]<-A[2]
	# Extend table to 100+ if not already 100+
	if(NFIN != 100 | nfnout != 100)  {
		oag[1:9]<-ARRAY[NEND,1:9]
		ictr<-6
		idx<-2
		oagL<-oag[5]
		for (i in seq(NEND,96,by=5))  {
			ictr<-ictr+1
			idx<-idx+5
			ARRAY[i,5]<-ARRAY[idx,5]
			oagL<-oagL-ARRAY[i,5]
			ARRAY[i,2]<-QVAL[ictr]
			ARRAY[i,3]<-ARRAY[i-5,3]*(1.0-ARRAY[i-5,2])
			ARRAY[i,4]<-ARRAY[i,3]*ARRAY[i,2]
			ARRAY[i,1]<-ARRAY[i,4]/ARRAY[i,5]
			ARRAY[i,9]<-G[ictr]
		}
		ARRAY[101,5]<-oagL
		ARRAY[101,2]<-1.0
		ARRAY[101,3]<-ARRAY[96,3]*(1.0-ARRAY[96,2])
		ARRAY[101,4]<-ARRAY[101,3]*ARRAY[101,2]
		ARRAY[101,1]<-ARRAY[101,4]/ARRAY[101,5]
		ARRAY[101,9]<-1.0/ARRAY[101,1]
		ARRAY[101,7]<-oagL
		ARRAY[101,8]<-ARRAY[101,7]/ARRAY[101,3]
		ARRAY[101,6]<-0.0
		idx<-101
		for (i in seq(NEND,96,by=5))  {
			idx<-idx-5
			ARRAY[idx,7]<-ARRAY[idx+5,7]+ARRAY[idx,5]
			ARRAY[idx,8]<-ARRAY[idx,7]/ARRAY[idx,3]
		}
		ARRAY[96,6]<-ARRAY[101,7]/ARRAY[96,7]
		nend5<-NEND-5
		for (i in seq(nend5,91,by=5)) {ARRAY[i,6]<-ARRAY[i+5,5]/ARRAY[i,5]}
		# New code allows output table to have flexible open age group, independent of input open age group (max of 100+).
		nend2<-nfnout+1
		ARRAY[nend2,2]<-1.0
		ARRAY[nend2,4]<-ARRAY[nend2,3]
		ARRAY[nend2,5]<-ARRAY[nend2,7]
		ARRAY[nend2,1]<-ARRAY[nend2,4]/ARRAY[nend2,5]
		ARRAY[nend2,9]<-1.0/ARRAY[nend2,1]
		ARRAY[nend2-5,6]<-ARRAY[nend2,7]/ARRAY[nend2-5,7]
	}  #  end of block to extend age group to 100+
	#  Might need output open age group (but can be determined from input parameters)
	#  Might need second output table because ARRAY calculates ages 0,1,2,3,4,5 then index every 5 years.
	
	
	
	
	#  for (ii in 1:101) {
	#    cat("\n",spaces(3),sprintf("%2i  %10.5f%10.5f%10.0f%10.0f%10.0f%10.5f%10.0f%10.3f%10.3f",ii-1,ARRAY[ii,1],ARRAY[ii,2],ARRAY[ii,3],ARRAY[ii,4],ARRAY[ii,5],ARRAY[ii,6],ARRAY[ii,7],ARRAY[ii,8],ARRAY[ii,9]))
	
	#  }
	
	
	
	
	return(list("ARRAY"=ARRAY,"ErrorCode"=ErrorCode))  
}


#' wrapper to a rote translation of the Abacus version LIFTB
#' @description A minimal argument wrapper that calls a minimally modified \code{AbacusLIFTB()}, an R transaltion of by UN staff of the Mortpak Fortran code. A wrapper is necessary because \code{AbacusLIFTB()} abides by some non-conventional indexing. We also take care of dimension naming.
#' @details The wrapper itself does little argument checking, but \code{AbacusLIFTB()} itself has nuanced argument checking and error reporting. Such errors rarely stop the function. Instead they are caught and returned. The wrapper returns such messages as warnings, but does not pass them on as objects. This function is used internally by \code{LTabr()}. Although this is a wrapper, it is not intended to be used directly.
#' @param Mx numeric vector of Mx in standard abridged age classes
#' @param qx numeric vector of qx in standard abridged age classes
#' @param mx_ind logical indicator of whether Mx or qx is given/preferred to be used
#' @param OAnew the desired open age group, max 100
#' @param Sex \code{"m"} or \code{"f"}. Anything not \code{"m"} treated as \code{"f"}
#' @export
#' @references 
#' \insertRef{mortpak1988}{DemoTools}
AbacusLIFTB_wrap <- function(Mx, qx, mx_ind = TRUE, OAnew = 100, Sex = "m"){
	
	if (missing(Mx) & !missing(qx)){
		mx_ind   <- FALSE
	}
	if (missing(qx) & !missing(Mx)){
		mx_ind   <- TRUE
	}
	if (mx_ind){
		mx_or_qx <- Mx
	} else {
		mx_or_qx <- qx
	}
	
	stopifnot(OAnew <= 100)
	N          <- length(mx_or_qx)
	AgeInt     <- inferAgeIntAbr(vec = mx_or_qx)
	OA         <- sum(AgeInt) - AgeInt[N]
	QXMX       <- rep(0,sum(AgeInt)+1)
	AgeI       <- cumsum(AgeInt)-AgeInt+1
	QXMX[AgeI] <- mx_or_qx
	NSEX       <- ifelse(Sex == "m",1,2)
	NTYPE      <- ifelse(mx_ind,2,1)
	OUT        <- AbacusLIFTB(NFIN = OA, NUMOUT = OAnew, 
			NTYPE = NTYPE, NSEX = 1, QXMX = QXMX)
	ARRAY           <- OUT$ARRAY
	AgeI            <- maxA2abridged(OAnew) + 1
	ARRAY           <- ARRAY[AgeI, ]
	if (sum(ARRAY == 0) / length(ARRAY) > .9){
		warning(OUT$ErrorCode)
	}
	colnames(ARRAY) <- c("Mx","qx","lx","dx","Lx","Sx","Tx","ex","ax")
	rownames(ARRAY) <- AgeI - 1
	ARRAY
}

