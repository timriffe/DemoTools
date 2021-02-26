# Author: ilya

# TR: TODO 1) make sure no negatives produced in value columns
# 2) also qx in range 0-1
# IK: done
###############################################################################
context("test-lt_abridged")


# testing function --------------------------------------------------------

lt_test_all_positive_plus_qx_lt_1 <- function(LT) {
    # check positive values
    expect_equal(
        LT %>% 
          # TR: open age AgeInt is NA, not handled well with logicals
            select(-c(Age,AgeInt)) %>% 
          # TR: rm is_weakly_less_than() since final 0 is valid sometimes.
            is_less_than(0) %>% #
            sum(),
        0
    )
    
    # check qx less than 1
    expect_equal(
        LT %>% 
            select(nqx) %>% 
            is_greater_than(1) %>% 
            sum(),
        0
    )
}


# PAS LTPOPDTH ------------------------------------------------------------

test_that("lt_abridged works on PAS example", {
    # trial code from PAS LTPOPDTH, North, Males, IMR = .1
    Exposures <- c(100958, 466275, 624134, 559559, 446736, 370653, 
                   301862, 249409, 247473, 223014, 172260, 149338, 
                   127242, 105715, 79614, 53660, 31021, 16805, 
                   8000, 4000, 2000, 1000)
    
    Deaths <- c(8674, 1592, 618, 411, 755, 1098, 1100, 1357,
                1335, 3257, 2200, 4023, 2167, 4578, 2956, 4212,
                2887, 2351, 1500, 900, 500, 300)
    # lower age bounds
    Age    <- c(0, 1, seq(5, 100, by = 5))
    AgeInt <- c(diff(Age), NA)
    
    PASLT <- lt_abridged(Deaths = Deaths,
                   Exposures = Exposures,
                   Age = Age,
                   AgeInt = AgeInt,
                   axmethod = "pas",
                   IMR = .1,
                   region = "n",
                   Sex = "m",
                   a0rule="cd")
    
    # The unsmoothed output from PAS
    # spreadsheet (allow for rounding error in last decimal)
    excheck <- c(56.31, 61.53, 58.35, 53.63,
                 48.81, 44.21, 39.83, 35.52,
                 31.43, 27.22, 24.09, 20.52,
                 18.12, 14.51, 12.42, 9.45, 7.85,
                 6.09, 4.95, 4.28, 3.85, 3.33)
    
    # final ex values are different because lifetable
    # closeout is via lifetable extention rather than mx inverse.
	# TR: closeout no longer comparable, only check first 19 diffs
    ind <- 1:19
    
    # test
    expect_equal(
        PASLT$ex[ind],
        excheck[ind],
        tolerance = .01
    )
    
    # positive, qx =< 1
    PASLT %>% lt_test_all_positive_plus_qx_lt_1()

})


# UN 1982 (p. 34) example -------------------------------------------------

test_that("lt_abridged works on UN 1982 (p. 34) example", {
    
    Mx <- c(.23669, .04672, .00982, .00511, .00697, .01036, .01169,
            .01332, .01528, .01757, .02092, .02517, .03225, .04241, 
            .06056, .08574, .11840, .16226, .23745)
    
    excheckUN <- c(35.000, 42.901, 47.190, 44.438, 40.523, 36.868, 
                   33.691, 30.567, 27.500, 24.485, 21.504, 18.599,
                   15.758, 13.080, 10.584, 8.466, 6.729, 5.312, 4.211)
    Age    <- c(0,1,seq(5,85,by=5))
    AgeInt <- age2int(Age, OAvalue = 5)
    
    # generate two variants: with and without PG's variants for ages 5-14
    UNLT1 <- lt_abridged(nMx = Mx,
                   Age = Age,
                   AgeInt = AgeInt,
                   axmethod = "un",
                   Sex = "m",
                   mod = FALSE)
    
    UNLT2 <- lt_abridged(nMx = Mx,
                   Age = Age,
                   AgeInt = AgeInt,
                   axmethod = "un",
                   Sex = "m",
                   mod = TRUE)
    
	# TR 2 Dec 2018: lt_abridged closeout is permanently
	# changed, influencing all ages, especially old ages.
	# results should be close in young ages
	ind <- 1:14
    # test
    expect_equal(
        UNLT1$ex[ind],
        excheckUN[ind],
        tolerance = .02
    )
    
    # positive, qx =< 1
    UNLT1 %>% lt_test_all_positive_plus_qx_lt_1()
    UNLT2 %>% lt_test_all_positive_plus_qx_lt_1()

})



# Mortpak example --------------------------------------------------------

# TR Oct 13: I've been tinkering with these, and for now we no longer
# have the option of reproducing these results. This is because we have 
# an internally consistent closeout procedure that does not include
# the ABACUS method as an option. It wins us consistency, though.
# IFF I revitalize the Mortpak closeout method then it'll come back
# in via args and be tested, but for now can't do that.
test_that("lt_abridged works on Mortpak example (United Nations 1988, p. 82)", {
    
    MPnMx <- c(0.12846, 0.02477, 0.00603, 0.0034,
               0.00417, 0.00513, 0.00581, 0.00645, 0.00725,
               0.00813, 0.00913, 0.01199, 0.01647,
               0.0256, 0.04047, 0.06624, 0.10638, 0.19611)
    
    MPexcheck <- c(49.997, 55.675, 57.245, 53.921,
                   49.803, 45.799, 41.922, 38.084, 34.249,
                   30.420, 26.578, 22.701, 18.945,
                   15.349, 12.095, 9.240, 6.903, 5.099)
    Age <-c(0, 1, seq(5, 80, by = 5))
    # First with lifetable extention to 100
    MP_UNLT100 <- lt_abridged(
        nMx = MPnMx,
        Age = Age,
        AgeInt = inferAgeIntAbr(vec = MPnMx),
        axmethod = "un",
        Sex = "f",
        mod = FALSE,
        OAnew = 100
    )
    
    # lifetable to original open age group
    MP_UNLT80 <- lt_abridged(
        nMx = MPnMx,
        Age = c(0, 1, seq(5, 80, by = 5)),
        AgeInt = age2int(Age),
        axmethod = "un",
        Sex = "f",
        mod = FALSE,
        OAnew = 80
    )
    
    # same, but truncated at 60
    MP_UNLT60 <- lt_abridged(
        nMx = MPnMx,
        Age = c(0, 1, seq(5, 80, by = 5)),
        AgeInt = inferAgeIntAbr(vec = MPnMx),
        axmethod = "un",
        Sex = "f",
        mod = FALSE,
        OAnew = 60
    )
    
    # tests
    
    # matches published results closely
# TR: Oct 13 removed this test, although keep an eye on it.
# the mortpak extrpolation has been removed and replaced with 
# MortalityLaws extrapolation, which is more flexible. Some things
# still under consideration, but for now this indeed shouldn't be
# the same (which is odd I know).
#    expect_equal(
#        MP_UNLT100$ex[1:length(MPexcheck)],
#        MPexcheck,
#        tolerance = .002
#    )
    
# TR Oct 13 2018: now removed this one as we no longer have a comparable case!
# new procedure is to extrapolate to a very high age ALWAYS using
# MotalityLaws, then truncate to OAnew, 
 #   expect_equal(
 #       MP_UNLT80$ex,
 #       MPexcheck,
 #       tolerance = .002
 #   )
    
    # identical results irrespective of max age
    expect_equal(
        MP_UNLT60$ex,
        MP_UNLT80$ex[1:length(MP_UNLT60$ex)],
        tolerance = 1e-12
    )
    
	# TR: as of now closeout is always the same and there
    # is no way to turn it off.
    expect_equal(
        MP_UNLT80$ex,
        MP_UNLT100$ex[1:length(MP_UNLT80$ex)],
        tolerance = 1e-12
    )
    
    # positive, qx =< 1
    MP_UNLT100 %>% lt_test_all_positive_plus_qx_lt_1()
    MP_UNLT80 %>% lt_test_all_positive_plus_qx_lt_1()
    MP_UNLT60 %>% lt_test_all_positive_plus_qx_lt_1()
})
