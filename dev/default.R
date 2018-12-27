
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

# agesmth recommended
h1 <- heapify(smoothed, Age, p0 = 1, p5 = .5, 20, 95)
lines(Age, h1)
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

# case where beers() sufficient, but sprague not enough:
h2 <- heapify(smoothed, Age, p0 = .5, p5 = .5, 20, 95)
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
h3 <- heapify(smoothed, Age, p0 = .2, p5 = .05, 20, 95)
plot(Age, h3, col = "red", type='l')
lines(Age, sprague(h3,Age,OAG=FALSE),lwd=2)
# beers not necessary
lines(Age, beers(h3,Age,OAG=FALSE),col="blue",lwd=2)

# very low roughness:
five_year_roughness(h3,Age)  # maaybe
zero_pref_sawtooth(h3,Age)   # none detected.
Myers(h3,Age)
Whipple(h3,Age)
Spoorenberg(h3,Age)
