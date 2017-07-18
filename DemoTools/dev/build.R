
# Author: tim
###############################################################################

library(devtools)
document("/home/tim/git/DemoTools/DemoTools")


load_all("/home/tim/git/DemoTools/DemoTools")

#.Fortran("ABBREV",D1=rep(5,25),D2=rep(5,26))
#
#install.packages("inline")
#body = "C     SUBROUTINE INTRP(LOE,NWRIT,N,YA,YB,YC,A,C,B)
#C-----------------------------------------------------------------------
#C----- PROGRAM NO. 0520
#C-----------------------------------------------------------------------
#C----- (Modified April 17, 1985 for MS-DOS)
#C-----------------------------------------------------------------------
#C----- THE INPUT ARGUMENTS TO THIS SUBROUTINE ARE LOE, NWRIT, N, YA, YB,
#C----- YC, A, AND C.
#C----- THE OUTPUT ARGUMENT FROM THIS SUBROUTINE IS B.
#C----- LOE  INDICATES EITHER LINEAR OR EXPONENTIAL INTERPOLATION.
#C----- NWRIT  PRINT INDICATOR.  PRINTS WHEN NONZERO.
#C----- N  IS THE NUMBER OF INTERPOLATIONS TO BE PERFORMED.
#C----- YA  REFERENCE POINT FOR VALUES IN A.
#C----- YB  POINT TO BE INTERPOLATED FOR.
#C----- YC  REFERENCE FOR VALUES IN C.
#C----- A  CONTAINS THE VALUES TO BE USED FOR THE INTERPOLATION REFERRING
#C-----       TO POINT YA.
#C----- C  CONTAINS THE VALUES TO BE USED FOR THE INTERPOLATION REFERRING
#C-----       TO POINT YC.
#C----- B  VALUES THAT WERE INTERPOLATED FOR BY THE SUBROUTINE.  THEY
#C-----       PERTAIN TO THE POINT YB.
#C-----------------------------------------------------------------------
#C-----------------------------------------------------------------------
#      DIMENSION A(99),C(99),B(99)
#*-----
#      write(*,1)
#    1 format(5x,'Running SUBROUTINE INTRP ...')
#*-----
#      NPRNT=6
#C
#C-----------------------------------------------------------------------
#C----- ERROR CHECK
#C-----------------------------------------------------------------------
#C
#      IF( YC - YA ) 70,60,70
#   60 WRITE(NPRNT,222)YA
#  222 FORMAT(/,1X,51H*** INTRP ERROR NO. 0521 -- INPUT ERROR IN YA OR YC
#     *,/,1X,49H*** YEARS INTERPOLATED BETWEEN ARE BOTH EQUAL TO ,F8.3,
#     *38H THEY MUST NOT BE EQUAL TO EACH OTHER.)
#   63 WRITE(NPRNT,11)
#   11 FORMAT(/,5X,18HINTRP INPUT VALUES,/)
#      DO 65 I =1,N
#   65 B(I) = 0.0
#      WRITE(NPRNT,22) LOE,NWRIT,N,YA,YB,YC
#   22 FORMAT(5X,6HLOE = ,I1,/,5X,8HNWRIT = ,I1,/,5X,4HN = ,I2,
#     */,5X,5HYA = ,F13.5,/,5X,5HYB = ,F13.5,/,5X,5HYC = ,F13.5)
#      WRITE(NPRNT,33) (A(I),I=1,N)
#   33 FORMAT(/,5X,4HA = ,5(F13.5,3X),/,9X,5(F13.5,3X))
#      WRITE(NPRNT,44) (C(I),I=1,N)
#   44 FORMAT(/,5X,4HC = ,5(F13.5,3X),/,9X,5(F13.5,3X))
#      GO TO 1000
#C
#C-----------------------------------------------------------------------
#C----- CALCULATE DISTANCE BETWEEN POINTS AND DETERMINE WHETHER LINEAR
#C----- INTERPOLATION IS TO BE DONE.
#C-----------------------------------------------------------------------
#C
#   70 Z = ( YB - YA ) / ( YC - YA )
#      IF( LOE - 2 ) 80,160,80
#C
#C-----------------------------------------------------------------------
#C----- LINEAR INTERPOLATION
#C-----------------------------------------------------------------------
#C
#   80 DO  90 I = 1,N
#   90 B(I) = A(I) - Z * ( A(I) - C(I) )
#  120 CONTINUE
#      IF( LOE - 1 ) 130,140,130
#  130 WRITE(NPRNT,444)LOE
#  444 FORMAT(/,1X,19H*** INTERP WARNING:,/,1X,25H*** INTERPOLATION CODE
#     * IS,I2,85H IT MUST BE 1 (LINEAR) OR 2 (EXPONENTIAL); BOTH LINEAR A
#     *ND EXPONENTIAL WERE COMPUTED.,
#     */,1X,52H*** EXPONENTIALLY INTERPOLATED VALUES WERE RETURNED.)
#  140 CONTINUE
#      IF( NWRIT ) 150,1000,150
#  150 WRITE(NPRNT,555)YA,YB,YC
#  555 FORMAT(/,12X,20HLINEAR INTERPOLATION,//,13X,18HINTERPOLATED VALUE,
#     *//,1X,3(F13.5,6X),/)
#      DO 155 I = 1,N
#  155 WRITE(NPRNT,777) A(I),B(I),C(I)
#      IF (LOE-1)160,1000,160
#C
#C-----------------------------------------------------------------------
#C----- EXPONENTIAL INTERPOLATION
#C-----------------------------------------------------------------------
#C
#  160 DO 170 I = 1,N
#      IF (ABS(A(I)) - 0.0000001) 163,163,170
#  163 WRITE(NPRNT,55)
#   55 FORMAT(/,1X,44H*** INTRP ERROR NO. 0522 -- INPUT ERROR IN A
#     *,/,1X,  '*** ONE VALUE OF A IS ZERO, EXPONENTIAL INTERPOLATION IS
#     *IMPOSSIBLE. ')
#      GO TO 63
#  170 B(I) = A(I) * EXP(ALOG( C(I) / A(I) ) * Z )
#      IF( NWRIT) 180,1000,180
#  180 WRITE(NPRNT,666)YA,YB,YC
#  666 FORMAT(/,13X,25HEXPONENTIAL INTERPOLATION,//,17X,18HINTERPOLATED V
#     *ALUE,//,1X,3(F13.5,6X),/)
#  190 DO 200 I = 1,N
#  200 WRITE(NPRNT,777)A(I),B(I),C(I)
#  777 FORMAT(1X,3(F13.5,6X))
# 1000 RETURN
#   "
#library(inline)
#cfunction(signature(LOE="integer",
#				NWRIT="integer",
#				N = "integer",
#				YA = "double",
#				YB = "double",
#				YC = "double",
#				A = "double",
#				C = "double",
#				B = "double"),body, convention=".Fortran")