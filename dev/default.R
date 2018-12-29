
# Author: tim
###############################################################################


# testing for default processing.

# take arbitrary data in single or other age groups, and produce single age data... somehow.

# 1) need opag(), with decent default standard.
library(DemoTools)
library(magrittr)
pop_in <- pop1m_pasex
Age    <- 0:99

smoothed <- pop_in %>%
		agesmth(Age=Age,
				method = "Strong", 
 				OAG = FALSE, 
 				young.tail = "Arriaga") %>%
        sprague(OAG=FALSE) %>%
		rescale.vector(scale=sum(pop1m_pasex))

# agesmth recommended, followed by sprague
h1 <- heapify(smoothed, Age, p0 = .1, p5 = .1, 20, 95)
# case where beers() sufficient, but sprague not enough:
h2 <- heapify(smoothed, Age, p0 = .2, p5 = .2, 20, 95)
# go even gentler: sprague-only sufficient.
h3 <- heapify(smoothed, Age, p0 = .3, p5 = .2, 20, 95)
# insane 0-pref, light 5-pref
h4 <- heapify(smoothed, Age, p0 = .5, p5 = .3, 20, 95)
# insane heaping, even 0-5 pref
h5 <- heapify(smoothed, Age, p0 = .43, p5 = .5, 20, 95)

cols <- RColorBrewer::brewer.pal(7,"Reds")
plot(Age, smoothed, type='l',lwd=4,col=gray(.7),ylab = "counts")
lines(Age, h1, col = cols[7])
lines(Age, h2, col = cols[6])
lines(Age, h3, col = cols[5])
lines(Age, h4, col = cols[4])
lines(Age, h5, col = cols[3])


assess_pop <- function(Pop, Age, ageMin = 30, ageMax = min(c(70,max(Age[Age%%10 == 0])-10))){
	
	# single age checks only if single age data given:
	wi <- NULL
	mi <- NULL
	do.single <- is_single(Age)
	if (do.single){
		wi <- Whipple(Value = Pop, Age = Age, ageMin = ageMin, ageMax = ageMax, digit = c(0,5))
		mi <- Myers(Value = Pop, Age = Age, ageMin = ageMin, ageMax = ageMax)
	}
	# in all cases can run 5-year checks:
	si <- as.numeric(zero_pref_sawtooth(Value = Pop, Age = Age, ageMin = ageMin, ageMax = ageMax))
	ri <- five_year_roughness(Value = Pop, Age = Age, ageMin = ageMin, ageMax = ageMax)
	
	# negative lower bound to capture 0s
	# whipple break points:
	w_breaks <- c(-.01, 1.05, 1.2, 1.5, 10)
# 4 levels: negligible | light | moderate | heavy
# Myers break points:
	m_breaks <- c(-.01, 1, 3, 10, 100)
# 4 levels: negligible | light | moderate | heavy
	
# five-year-roughness breaks:
	r_breaks <- c(-.01, .1, .5, Inf)
# 3-level: light | mod | heavy (if not met with some sawtooth, then needs visual confirmation.)
	
# 0-pref-sawtooth beaks
	s_breaks <- c(-.01,1, 1.1, 1.5, Inf) # last interval most aggressive
# 4 levels: negligible | light | mod | heavy
	
	r_level <- as.character(cut(ri, r_breaks, labels = c("light", "moderate", "heavy")))
	s_level <- as.character(cut(si, s_breaks, labels = c("negligible", "light", "moderate", "heavy")))
	
	if (do.single){
		w_level <- as.character(cut(wi, w_breaks, labels = c("negligible", "light", "moderate", "heavy")))
		m_level <- as.character(cut(mi, m_breaks, labels = c("negligible", "light", "moderate", "heavy")))
	} else {
		w_level <- NA
		m_level <- NA
	}

	c(Whipple = w_level, Myers = m_level, Roughness = r_level, Sawtooth = s_level)
}

# assessment coming from assess_pop()
plan_pop_adjustment <- function(assessment){
	
	levs <- c("negligible","light","moderate","heavy")
	need_plan <- TRUE
	# if we have single age data we either do or do not look 
	# at the 5-year data, depending on how it is
	if (!any(is.na(assessment[c("Whipple","Myers")]))){
		if (any(assessment[c("Whipple","Myers")] == "negligible")){
			the_plan  <- "as-is"
			need_plan <- FALSE
		}
		if (any(assessment[c("Whipple","Myers")] %in% c("light","moderate")) &
				assessment["Sawtooth"] == "negligible" &
				assessment["Roughness"] == "light"){
			the_plan  <- "sprague"
			need_plan <- FALSE
		}
		# can have terrible heaping with little consequence in 5-year age groups too
		if (any(assessment[c("Whipple","Myers")] == "heavy") &
				assessment["Sawtooth"] == "light"){
			the_plan  <- "agesmth(constrained)-sprague"
			need_plan <- FALSE
		}
	}
	
	# otherwise, the rest depends on the 5-year data assessment.
	
	# beers is the lightest treatment
	if (assessment["Roughness"] %in% c("light","moderate") & 
			assessment["Sawtooth"] == "light" & need_plan){
		    the_plan <- "beers"
			need_plan <- FALSE
	}
	if (assessment["Roughness"] %in% c("light","moderate") & 
			assessment["Sawtooth"] == "moderate" & need_plan){
		the_plan <- "agesmth(constrained)-sprague"
		need_plan <- FALSE
	}
	if ((assessment["Sawtooth"] == "heavy" | 
					(sum(assessment == "heavy") > 2)) & need_plan){
		the_plan <- "agesmth(aggressive)-sprague"
		need_plan <- FALSE
	}
	

	# for example if roughness is heavy, but sawtooth not triggered, what's going on?
	# maybe just abrupt fertiltiy history, but worth looking at.
	if (need_plan){
		the_plan <- "visual assessment required"
	}
	
	the_plan	
}

adjust_pop_with_plan <- function(Pop, Age, plan, OAG = TRUE){
	
	if (plan == "as-is"){
		return(Pop)
	}
	if (plan == "sprague"){
		out <- sprague(Pop, Age, OAG = OAG)
		if (any(out < 0)){
			cat("Warning: sprague produced negative values")
		}
		return(out)
	}
	if (plan == "beers"){
		out <- beers(Pop, Age, OAG = OAG)
		if (any(out < 0)){
			cat("Warning: beers produced negative values")
		}
		return(out)
	}

	if (plan == "agesmth(constrained)-sprague"){
#		assess_in <- assess_pop(Pop, Age)
		
		out1 <- agesmth(Pop, Age, method = "Carrier-Farrag")
		out  <- sprague(out1, Age = names2age(out1), OAG = OAG)
		if (any(out < 0)){
			cat("Warning: sprague produced negative values")
		}
		return(out)
#		out2 <- agesmth(Pop, Age, method = "KKN")
#		out3 <- agesmth(Pop, Age, method = "Zigzag")
#		
#		assess_out1 <- assess_pop(out1, names2age(out1))
#		assess_out2 <- assess_pop(out2, names2age(out2))
#		assess_out3 <- assess_pop(out3, names2age(out3))
#		
#		plot(names2age(out1), out1,type='o')
#		lines(names2age(out1), out2,col="red")
#		lines(names2age(out1), out3,col="blue")
		# need ordinal shortcut... here zigzag does better?
		# not really if you look at visual...
	}
	if (plan == "agesmth(aggressive)-sprague"){
		out1 <- agesmth(Pop, Age, method = "Strong")
		
		# this will affect tails, even if tails not adjusted!
		out2 <- rescale.vector(out1,scale = sum(Pop))
		out  <- sprague(out1, Age = names2age(out2), OAG = OAG)
		if (any(out < 0)){
			cat("Warning: sprague produced negative values")
		}
		return(out)
}
}

adjust_pop <- function(Pop, Age, OAG = TRUE, 
		ageMin = 30, ageMax = min(c(70, max(Age[Age%%10 == 0]) - 10))){
	assessment <- assess_pop(Pop, Age, ageMin = ageMin, ageMax = ageMax)
	plan       <- plan_pop_adjustment(assessment)
	if (plan == "visual assessment required"){
		cat(plan,"\n")
	}
	adjust_pop_with_plan(Pop = Pop, Age = Age, plan = plan, OAG = OAG)
}


plot(Age,h1,pch=1,cex=.7)
lines(Age,adjust_pop(h1,Age,FALSE),col="red")

plot(Age,h2,pch=1,cex=.7)
lines(Age,adjust_pop(h2,Age,FALSE),col="red")

plot(Age,h3,pch=1,cex=.7)
lines(Age,adjust_pop(h3,Age,FALSE),col="red")

plot(Age,h4,pch=1,cex=.7)
lines(Age,adjust_pop(h4,Age,FALSE),col="red")

plot(Age,h5,pch=1,cex=.7)
lines(Age,adjust_pop(h5,Age,FALSE),col="red")

plot(Age,pop1m_pasex,pch=1,cex=.7)
lines(Age,adjust_pop(pop1m_pasex,Age,FALSE),col="red")

plot(0:100,pop1m_ind,pch=1,cex=.7)
lines(0:100,adjust_pop(pop1m_ind,0:100,TRUE),col="red")

# ----------------------------------------------------
# below various pieces of scratch code used in testing.



# if single age data don't trigger Myers or Whipple, then no need for 5-year checks.
# If single age tripped, and roughness light and sawtooth negligible, then sprague enough
# if "                 , and roughness moderate and sawtooth light then beers enough
# if "                 , and roughness moderate and sawtooth mod then agesmth constrained
# if "                 , and roughness >mod | sawtooth > mod then agesmth aggressive

inda <- assess_pop(pop1m_ind,0:100)
plan_pop_adjustment(inda)

plot(0:100,pop1m_ind,type="o")
lines(0:100, sprague(agesmth(pop1m_ind,0:100,method="Strong"),Age=seq(0,100,by=5),OAG=TRUE))
# hmmmm. It would seem natural variation could easily bring Myers of .5 or Whipple of 1.03
# maybe best not to do anything.
# per UN Whipple < 1.05 very accurate. (Myers < .5)

# Sprague or nothing if:
(Whipple(h1,Age) < 1.05 | Myers(h1,Age) < 1)  &
( five_year_roughness(h1,Age) < .1 | zero_pref_sawtooth(h1,Age) == 0 )
h1a <- assess_pop(h1,Age)
plan_pop_adjustment(h1a)

# Beers if: not sprague... AND:
((Whipple(h2, Age) < 1.2 & Whipple(h2, Age) >= 1.05) | 
			(Myers(h2, Age) >= 1 & Myers(h2, Age) < 3)) &
		( five_year_roughness(h2, Age) < .1 | zero_pref_sawtooth(h2,Age) < 1.1)
h2a <- assess_pop(h2,Age)
plan_pop_adjustment(h2a)

# still didnt trigger beers
h2.2 <- heapify(smoothed, Age, p0 = .21, p5 = .2, 20, 95)
h2.2a <- assess_pop(h2.2,Age)
plan_pop_adjustment(h2.2a)

# agesmth(constrained)-sprague if:
((Whipple(h3,Age) >= 1.2 & Whipple(h3,Age) < 1.5) | 
			(Myers(h3,Age) >= 3) & (Myers(h3,Age) < 10) & # same
		( five_year_roughness(h3,Age) >= .1 & five_year_roughness(h3,Age) < .5 | 
			(zero_pref_sawtooth(h3,Age) > 1.05 & zero_pref_sawtooth(h3,Age) < 1.5 ))
h3a <- assess_pop(h3,Age)
plan_pop_adjustment(h3a)

# agesmth(strongly)-sprague if:
(Whipple(h4,Age) >= 1.5 | Myers(h4,Age) >= 10) &
		( five_year_roughness(h4,Age) >= .5 | zero_pref_sawtooth(h4,Age) >= 1.5 )
h4a <- assess_pop(h4,Age)
plan_pop_adjustment(h4a)

(Whipple(h5,Age) >= 1.5 | Myers(h5,Age) >= 10) & # same
		( five_year_roughness(h5,Age) >= .5 | zero_pref_sawtooth(h5,Age) < 1.5 )
h5a <- assess_pop(h5,Age)
plan_pop_adjustment(h5a)

plot(Age,smoothed,type='l')
lines(Age,beers(h5,Age,OAG=FALSE),col="blue")
lines(Age,sprague(agesmth(h5,Age,method="Carrier-Farrag"),Age=seq(0,95,by=5),OAG=FALSE),col="red")

# which ones constrained?
A5 <- seq(0,95,by=5)
# CF is 10-year constrained- no tails
groupAges(agesmth(h1,Age,method="Carrier-Farrag"),Age=A5,N=10) - groupAges(h1,Age,N=10)
# KKN is 10-year constrained- no tails
groupAges(agesmth(h1,Age,method="KKN"),Age=A5,N=10) - groupAges(h1,Age,N=10)
# Arriaga is 10-year constrained-0 includes tails
groupAges(agesmth(h1,Age,method="Arriaga"),Age=A5,N=10) - groupAges(h1,Age,N=10)
# UN is not constrained, even in total
sum(agesmth(h1,Age,method="United Nations")) - sum(h1)
groupAges(agesmth(h1,Age,method="United Nations"),Age=A5,N=5) - groupAges(h1,Age,N=5)
# Strong is not constrained, except in the 10-year young tails. old tail unchanged.
# penultimate 10-year age group is constrained.
sum(agesmth(h1,Age,method="Strong")) - sum(h1)
agesmth(h1,Age,method="Strong") - groupAges(h1,Age,N=5)
groupAges(agesmth(h1,Age,method="Strong"),Age=A5,N=10) - groupAges(h1,Age,N=10)
# Zigzag is constrained in tails (longer old tail), and total is constrained,
# mid middle ages not in constrained groups.
sum(agesmth(h1,Age,method="Zigzag")) - sum(h1)
agesmth(h1,Age,method="Zigzag") - groupAges(h1,Age,N=5)
groupAges(agesmth(h1,Age,method="Zigzag"),Age=A5,N=10) - groupAges(h1,Age,N=10)
plot(A5,agesmth(h4,Age,method="Zigzag",ageMax=90) )

sum(agesmth(h1,Age,method="mav")) - sum(h1)



# Beers if:
# Whipple
DemoTools::Whipple
A5 <- seq(0,95,by=5)
h1_5 <- groupAges(h1,Age=Age)
plot(A5, h1_5/5, type="s")
lines(Age, sprague(h1_5,A5,OAG=FALSE),col="blue", lwd = 2)
lines(Age, beers(h1_5,A5,OAG=FALSE),col="red", lwd = 2)
# still slightly wavy, but conclude beers is stronger than sprague.
plot(Age, smoothed, type='l')
lines(Age, beers(h1_5,A5,OAG=FALSE),col="red", lwd = 2)

# if we were going to merely group-graduate then beers preferred, (not constrained!)
groupAges(beers(h1_5,A5,OAG=FALSE),Age) - h1_5
# constrained!
groupAges(sprague(h1_5,A5,OAG=FALSE),Age) - h1_5

# -----------------------------------------------------------------------------------#
# in principle, if we use agesmth(), then always follow with sprague, not beers!     #
# i.e. keep age-group shifting all together in one step, not parsed into two steps.  #
# -----------------------------------------------------------------------------------#

five_year_roughness(h1,Age)  # probable yes
zero_pref_sawtooth(h1,Age)   # definite yes

h2_5 <- groupAges(h2, Age)

plot(A5, h2_5/5, type='s')
lines(A5+2.5, h2_5/5)
lines(Age, h2,col="red")

five_year_roughness(h2,Age)  # maaybe
zero_pref_sawtooth(h2,Age)   # none detected.
# conclude likely sprague or at most beers sufficient.

plot(Age, h2, col = "red", type='l')
lines(Age, sprague(h2,Age,OAG=FALSE),lwd=2,lty=2)
lines(Age, beers(h2,Age,OAG=FALSE),col="blue",lwd=2)

# -----------------------------
# go even gentler
plot(Age, h3, col = "red", type='l')
lines(Age, sprague(h3,Age,OAG=FALSE),lwd=2)
# beers not necessary
lines(Age, beers(h3,Age,OAG=FALSE),col="blue",lwd=2)





Ages         <- seq(0, 80, by = 5)

# names a bit flexible:
cf <- agesmth(Value = pop5m_pasex, 
		Age = Ages, 
		method = "Carrier-Farrag", 
		OAG = TRUE)
# old.tail be default same as young.tail
# "cf" also works

# no need to specify tails for Arriaga or Strong
arr <- agesmth(Value = pop5m_pasex, 
		Age = Ages, 
		method = "Arriaga", 
		OAG = TRUE)
strong <- agesmth(Value = pop5m_pasex, 
		Age = Ages, 
		method = "Strong", 
		OAG = TRUE)
# other methods:
un <- agesmth(Value = pop5m_pasex, 
		Age = Ages, 
		method = "United Nations", 
		OAG = TRUE)
kkn <- agesmth(Value = pop5m_pasex, 
		Age = Ages, 
		method = "Karup-King-Newton", 
		OAG = TRUE)
# zigzag, not plotted.
zz <- agesmth(pop5m_pasex,Ages,OAG=TRUE,method="Zigzag",ageMin = 30, ageMax = 80)
# mav, not plotted.
ma3 <- agesmth(pop5m_pasex,Ages,OAG=TRUE,method="MAV",n=3)

## Not run:

plot(Ages,pop5m_pasex,pch=16)
lines(Ages, cf)
lines(Ages, arr, col = "red")
lines(Ages, strong, col = "#FF000080", lwd = 3)
lines(Ages, kkn, col = "blue")
lines(Ages, un, col = "magenta")
legend("topright",
		pch=c(16,NA,NA,NA,NA,NA),
		lty = c(NA,1,1,1,1,1),
		lwd = c(NA,1,1,3,1,1),
		col = c("black","black","red","#FF000080","blue","magenta"),
		legend = c("orig 5","Carrier-Farrag","Arriaga","Strong","KKN","UN"))

# Strong > UN > Arriaga > CarrierFarrag

five_year_roughness(pop5m_pasex,Ages) # 5
zero_pref_sawtooth(pop5m_pasex,Ages)  # 1.6



 Myers(rep(c(rep(0,9),1),10),Age=0:99)
 Myers(rep(1,100),Age=0:99)
 Bachi(rep(1,100), Age=0:99, ageMin = 10, ageMax = 89)
 Bachi(rep(c(rep(0,9),1),10), Age=0:99, ageMin = 10, ageMax = 90)
 Noumbissi(rep(1,100),0:99,ageMin=20,ageMax=80,digit=0)
 Noumbissi(rep(c(rep(0,9),1),10),0:99,ageMin=20,ageMax=80,digit=9)
 Spoorenberg(rep(c(rep(0,9),1),10),0:99,ageMin=20,ageMax=80)
 Spoorenberg(rep(1,100),0:99,ageMin=20,ageMax=80)
 Whipple(rep(c(rep(0,9),1),10), Age=0:99, ageMin = 20, ageMax = 70,digit=9)
CoaleLi(pop1m_ind, 0:100, ageMin = 60, ageMax = max(Age), terms = 5, digit = 0)
CoaleLi(rep(1,100), 0:99, ageMin = 60, ageMax = max(Age), terms = 5, digit = 0)
CoaleLi(rep(c(rep(0,9),1),10), 0:99, ageMin = 60, ageMax = max(Age), terms = 5, digit = 5)
KannistoHeap(Value = pop1m_pasex, Age = 0:99, Agei = 90)
KannistoHeap(Value = pop1m_ind, Age = 0:100, Agei = 90)

 
 Myers(pop1m_pasex,0:99,ageMin=20,ageMax=80)
 Bachi(pop1m_pasex,0:99,ageMin=20,ageMax=80)
 
 Noumbissi(pop1m_pasex,0:99,ageMin=20,ageMax=80,digit=0)
 Whipple(pop1m_pasex,0:99,ageMin=20,ageMax=80,digit=0)
 CoaleLi(pop1m_ind, Age, ageMin = 60, ageMax = max(Age), terms = 5, digit = 0)
 P5 <- groupAges(pop1m_ind,0:100)
 A5 <- seq(0,100,by=5)
 zero_pref_sawtooth(pop1m_ind, Age = 0:100, ageMin = 40, ageMax = 90)
 
 five_year_roughness(agesmth(pop1m_ind, Age = 0:100,method="Strong"), Age = A5, ageMin = 40, ageMax = 80)
 
 