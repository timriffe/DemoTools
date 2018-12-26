
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


heapify <- function(Value, Age, p0=2, p5=p0, ageMin = 25, ageMax = max(Age[Age %% 5 == 0])){
	# pascal weights
	# x3,x4,x5,x6,x7
	fivepdf <- rescale.vector(c(1,4,6,4,1)) 
	#x7,x8,x9,x0,x1,x2,x3 
	#tenpdf  <- rescale.vector(c(1,6,15,20,15,6,1)) * p0
	
	# center ages:
	
	A10 <- Age[Age %% 10 == 0 & Age >= ageMin & Age <= ageMax]
	A5  <- Age[Age %% 5 == 0 & !Age%in%A10 & Age >= ageMin & Age <= ageMax]
	
	# rest could happen in a matrix, but this might be clearer to read:
	ai <- sort(c(A5,A10))
	for (a in ai){
		if (a %% 10 == 0){
			dist <- 2
			pdf.i <- fivepdf * p0
		} else {
			dist <- 2
			pdf.i <- fivepdf * p5
		}
		sendi <- Age >= (a - dist) & Age <= (a + dist)
		Vsend <- Value[sendi] * pdf.i
		Value[sendi] <- Value[sendi] - Vsend
		Value[Age == a] <- Value[Age == a] + sum(Vsend)
	}
	Value
}


cbind(pop1m_pasex,heapify(Value,Age,2,1.5,ageMin=20))


plot(Age,Value,type='l')
points(Age,heapify(Value,Age,2.5,1,ageMin=20))
points(Age,pop1m_pasex,pch='x')

VH  <- heapify(Value,Age,2,1.5,ageMin=20)
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

zero_pref_sawtooth <- function(Value, Age, ageMin = 30, ageMax = max(Age[Age %% 10 == 0])){
	VH5    <- groupAges(Value,Age,5)
	A5     <- names2age(VH5)
	adj2   <- avg_adj(VH5)
	ai     <- A5 >= ageMin & A5 <= ageMax
	m05    <- suppressWarnings(matrix((VH5 / adj2)[ai], nrow = 2))
	m05    <- m05[, -ncol(m05)]
	
	# need rather consistent x0 > x5 pattern
	if (sum(diff(sign(log(m05))) == -2) < (ncol(m05)-2)){
		return(FALSE)
	}
	
	1/ratx(rowMeans(m05,na.rm=TRUE))
}
zero_pref_sawtooth

zero_pref_sawtooth(Value,Age)
zero_pref_sawtooth(pop1m_pasex,Age,ageMin=40,ageMax=60)
zero_pref_sawtooth(heapify(Value,Age,2,1.5,ageMin=20),Age)

Value <- heapify(Value,Age,3,1.5,ageMin=20)

#matplot(0:9,matrix(exp(abs(log(pop1m_pasex / Value))),nrow=10,byrow=TRUE),type='l')