
# Author: tim
###############################################################################
context("IRDID")

# Table 7.20 in Siegel and Swanson (2004). Make sure cited in IRDID.R
# under @references tag (using \insertRef{bibtextag}{DemoTools}

test_that("ID and IRD reproduce source examples", {
			pop1 <- c(7.38,14.16,14.79,17.36,15.11,10.14,8.50,7.28,5.28)
			pop2 <- c(6.48,12.27,15.25,15.10,14.66,10.80,8.95,9.28,7.21)
			
			# this is a set of expectations that are checked. The numbers 
			# 5.5 and 6.7 are pulled from the book example.
			expect_equal(ID(pop1, pop2), 5.5, tolerance = 0.01)
			expect_equal(IRD(pop1, pop2), 6.7, tolerance = 0.01)
			}
)

test_that("ID and IRD scale invariant", {
			# TR: it seems redundant to have to re-assign this, but the
			# tests appeared to break if I didn't. 
			pop1 <- c(7.38,14.16,14.79,17.36,15.11,10.14,8.50,7.28,5.28)
			pop2 <- c(6.48,12.27,15.25,15.10,14.66,10.80,8.95,9.28,7.21)
			
			expect_equal(ID(pop1, pop2), ID(pop1, pop2*2))
			expect_equal(IRD(pop1, pop2), IRD(pop1, pop2*2))
		}
)

test_that("ID and IRD identical", {
			pop1 <- c(7.38,14.16,14.79,17.36,15.11,10.14,8.50,7.28,5.28)
			pop2 <- c(6.48,12.27,15.25,15.10,14.66,10.80,8.95,9.28,7.21)
			
			expect_equal(ID(pop1, pop1),0)
			expect_equal(IRD(pop1, pop1),0)
		}
)
