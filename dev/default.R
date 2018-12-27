
# Author: tim
###############################################################################


# testing for default processing.

# take arbitrary data in single or other age groups, and produce single age data... somehow.

# 1) need opag(), with decent default standard.
library(magrittr)
pop_in <- pop1m_pasex
Age <- 0:99

smoothed <- pop_in %>%
		agesmth(Age=Age,
				method = "Strong", 
 				OAG = FALSE, 
 				young.tail = "Arriaga") %>%
        sprague(OAG=FALSE) %>%
		rescale.vector(scale=sum(pop1m_pasex))
    
plot(Age, smoothed, type='l')

# agesmth recommended, followed by sprague
h1 <- heapify(smoothed, Age, p0 = 1, p5 = .5, 20, 95)
# case where beers() sufficient, but sprague not enough:
h2 <- heapify(smoothed, Age, p0 = .5, p5 = .5, 20, 95)
# go even gentler: sprague-only sufficient.
h3 <- heapify(smoothed, Age, p0 = .2, p5 = .05, 20, 95)
# insance 0-pref, light 5-pref
h4 <- heapify(smoothed, Age, p0 = 1.5, p5 = .5, 20, 95)

# agesmth recommended
plot(Age, smoothed, type='l',col=gray(.8),lwd=4)
lines(Age, h1, col = gray(.5))
lines(Age, h2, col = gray(.3))
lines(Age, h3, col = gray(.1))
lines(Age, h4, col = gray(.1))
# what does Myers say?
Myers(smoothed,Age) # .3 % would need to be redistributed
Myers(h1,Age) # 9 % would need to be redistributed
Myers(h2,Age) # 4% would need to be redistributed
Myers(h3,Age) # 1/2 % would need to be redistributed
DemoTools::Myers
# ratio (exaggeration)
Whipple(h1,Age) # 48 %
Whipple(h2,Age) # 25 %
Whipple(h3,Age) # 2.8 %

# hmmmm. It would seem natural variation could easily bring Myers of .5 or Whipple of 1.03
# maybe best not to do anything.
# per UN Whipple < 1.05 very accurate. (Myers < .5)

# Sprague or nothing if:
(Whipple(h3,Age) < 1.05 | Myers(h3,Age) < 1)  &
( five_year_roughness(h3,Age) < .1 | zero_pref_sawtooth(h3,Age) == 0 )

# Beers if: not sprague... AND:
(Whipple(h2,Age) < 1.25 | Myers(h2,Age) < 5)  &
		( five_year_roughness(h2,Age) < .2 | zero_pref_sawtooth(h2,Age) < .05)

# agesmth(constrained)-sprague if:
(Whipple(h1,Age) >= 1.25 | Myers(h1,Age) >= 5) &
		( five_year_roughness(h1,Age) >= .2 | zero_pref_sawtooth(h1,Age) > 0 )

# agesmth(strongly)-sprague if:
(Whipple(h4,Age) >= 1.25 | Myers(h4,Age) >= 5) &
		( five_year_roughness(h4,Age) >= .5 | zero_pref_sawtooth(h4,Age) > 1 )


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