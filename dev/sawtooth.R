
# Author: tim
###############################################################################

# similar to Feeney zigzag, need nice way to detect something of the form
# heaping worse on 10s than on 5s, ergo age groups X0 > X5. Both?
# This can be a yes-no step for whether to do smoothing.

Value <- pop1m_pasex
Age <- 0:99

plot(Age, Value)

# avg ratio?
V5 <- groupAges(Value, Age, N = 5)
A5 <- names2age(V5)
plot(A5,V5)

X  <- V5
XX <- matrix(X,ncol=2,byrow=TRUE)
shift.vector(XX[,2],-1)


age10 <- seq(20, 70, by = 10)
WI10  <- sapply(age10, Jdanov, Value = pop1m_pasex, Ages = Age)


n <- length(x)
avg_adj_pow <- function(x, pow = 2){
	(shift.vector(x, -1, NA) ^ pow + 
	 shift.vector(x, 1, NA) ^ pow) ^ (1 / pow)
}
V5 / avg_adj_pow(V5,10)

XX[, 1] / ma(XX[,2],2)
ma(XX[,1],2) / ma(XX[,2],2)

avg2  <- avg_adj(V5) 
V5 / avg2

# difference squared
diffsq <- (Smoothed - avg2)^2 * id

age10 <- seq(20, 70, by = 10)
WI10  <- sapply(age10, Jdanov, Value = pop1m_pasex, Age = 0:99)

Jdanov(Value=pop1m_pasex,Age=0:99,Agei=50)

age5 <- seq(25, 75, by = 10)
WI5  <- sapply(age5, Jdanov, Value = pop1m_pasex, Age = 0:99)

plot(age10, WI10, type = "o", main = "heaping on 0s gets worse over age in this census")
lines(age5,WI5)
#       1
#      1 1
#     1 2 1
#    1 3 3 1
#   1 4 6 4 1        # 5s from this line
# 1 5 10 10 5 1
#1 6 15 20 15 6 1    # 10s from this line

plot(0:4,Value[Age >=60 & Age < 65])
# by this method, heaping is a force of sending, not receiving
sum(Value)
Value <- sprague(agesmth(pop1m_pasex,Age,method = "Strong",OAG=FALSE,young.tail="Arriaga"),OAG=FALSE)



cbind(pop1m_pasex,heapify(Value,Age,2,1.5,ageMin=20))

Age <- 0:99
plot(Age,Value,type='l')
points(Age,heapify(Value,Age=0:99,1.8,1.1,ageMin=20))
points(Age,pop1m_pasex,pch='x')

VH  <- heapify(Value,Age,2,1.5,ageMin=20)
sum(Value)

VH5 <- groupAges(VH,Age,5)

plot(names2age(VH5),VH5)

# This is a sawtoooth situation
A5     <- names2age(VH5)
adj2   <- avg_adj(VH5)

# 0 preference > 5 preference
mean(exp(abs(log(VH5 / adj2))[A5 >= ageMin & A5 <= ageMax]),na.rm=TRUE)

#zero_pref <- function(Value, Age, ageMin = 25, ageMax = max(Age[Age %% 5 == 0])){
#	VH5    <- groupAges(Value,Age,5)
#	A5     <- names2age(VH5)
#	adj2   <- avg_adj(VH5)
#	ai     <- A5 >= ageMin & A5 <= ageMax
#	mean(exp(abs(log(VH5 / adj2))[ai]),na.rm=TRUE)
#}

#' detect if heaping is worse on terminal digit 0s than on 5s
#' @description Ages ending in 0 (0s) often have higher apparent heaping than ages ending in 5 (5s). If heaping ocurrs in roughly the same amount on 0s and 5s, then it may be sufficient to group data into 5-year age groups and then graduate back to single ages. However, if heaping is worse on 0s, then this procedure tends to produce a wavy pattern in count data, with 10-year periodicity. In this case it is recommended to use one of the methods of \code{agesmth()} as an intermediate step before graduation. 
#' @details Data is grouped to 5-year age bins. The ratio of each value to the average of its neighboring values is calculated. If 0s have stronger attraction than 5s then we expect these ratios to be >1 for 0s and <1 for 5s. Ratios are compared within each 10-year age group in the evaluated age range. If in the evaluated range there are at most two exceptions to this rule (0s>5s), then the ratio of the mean of these ratios is returned, and it is recommended to use a smoother method. Higher values suggest use of a more aggressive method. This approach is only slightly different from that of Feeney, as implemented in the \code{zigzag()} functions. This is not a general measure of roughness, but rather an indicator of this particular pattern of age attraction. 
#' @export
#' @inheritParams DemoTools heapify
#' @param ageMin integer evenly divisible by 10. Lower bound of evaluated age range, default 40.
#' @param ageMax integer evently divisibly by 5. Upper bound of evaluated age range, defaults to highest age evenly divisible by 10.
#' @return \code{FALSE} if sawtooth pattern is not detected, numeric otherwise.
#' @references
#' \insertRef{feeney1979}{DemoTools}
#' @examples
#' smoothed <- sprague(
#' 		agesmth(pop1m_pasex, 
#' 				Age, 
#' 				method = "Strong", 
#' 				OAG = FALSE, 
#' 				young.tail = "Arriaga"),
#' 		OAG = FALSE)
#' # not saw-tooth jagged
#' zero_pref_sawtooth(smoothed, Age)
#' # saw-tooth pattern detected in older ages
#' zero_pref_sawtooth(pop1m_pasex, Age)
#' # heaped, but no 0>5 preference
#' h1 <- heapify(smoothed, Age, p0 = 1, p5 = 1)
#' # heaping progressively worse on 0s than on 5s.
#' h2 <- heapify(smoothed, Age, p0 = 1.2, p5 = 1)
#' h3 <- heapify(smoothed, Age, p0 = 1.5, p5 = .8)
#' h4 <- heapify(smoothed, Age, p0 = 2, p5 = .5)
#' 
#' \dontrun{
#' 	plot(Age, smoothed, type='l')
#' 	lines(Age, h1,col="blue")
#' 	lines(Age, h2,col="green")
#' 	lines(Age, h3,col="red")
#' }
#' zero_pref_sawtooth(h1, Age)
#' # 0-preference > 5 pref
#' zero_pref_sawtooth(h2, Age)
#' # increasing values
#' zero_pref_sawtooth(h3, Age)
#' zero_pref_sawtooth(h4, Age,ageMin=35)

zero_pref_sawtooth <- function(Value, Age, ageMin = 40, ageMax = max(Age[Age %% 5 == 0])){
	
	# rather than stopifnot() check, just make it work.
	ageMin <- ageMin - ageMin %% 10
	# group to 5-year ages if not already.
	VH5    <- groupAges(Value, Age, 5)
	A5     <- names2age(VH5)
	# elementwise matched to avg of adjacent values
	adj2   <- avg_adj(VH5)
	# evalauated age range
	ai     <- A5 >= ageMin & A5 <= ageMax
	# matrix of ratios, 0s in row 1, 5s in row 2
	m05    <- suppressWarnings(matrix((VH5 / adj2)[ai], nrow = 2))
	# ensure no recycled values used
	if (sum(ai) %% 2 != 0){
		m05    <- m05[, -ncol(m05)]	
	}
	
	# need rather consistent x0 > x5 pattern
	if (sum(diff(sign(log(m05))) == -2) < (ncol(m05) - 2)){
		# i.e. it could be very rough still, but not necessarily
		# a regular sawtooth pattern. Possibly a visual assessment needed,
		# since could be real pattern, or other pattern of roughness that also
		# requires smoothing.
		return(FALSE)
	}
	
	# mean of 0s divided by mean of 5s, that simple.
	1 / ratx(rowMeans(m05, na.rm = TRUE))
}


