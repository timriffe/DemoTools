# Author: ilya
###############################################################################
context("test-ltpopdth")


# PAS LTPOPDTH ------------------------------------------------------------

test_that("LTabr works on PAS example", {
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
    
    PASLT <- LTabr(Deaths = Deaths,
                   Exposures = Exposures,
                   Age = Age,
                   AgeInt = AgeInt,
                   axmethod = "pas",
                   IMR = .1,
                   region = "n",
                   Sex = "m")
    
    # The unsmoothed output from PAS
    # spreadsheet (allow for rounding error in last decimal)
    excheck <- c(56.31, 61.53, 58.35, 53.63,
                 48.81, 44.21, 39.83, 35.52,
                 31.43, 27.22, 24.09, 20.52,
                 18.12, 14.51, 12.42, 9.45, 7.85,
                 6.09, 4.95, 4.28, 3.85, 3.33)
    
    # final ex values are different because lifetable
    # closeout is via lifetable extention rather than mx inverse.
    ind <- 1:20
    
    # test
    expect_equal(
        PASLT$ex[ind],
        excheck[ind],
        tolerance = .01
    )
})


# UN 1982 (p. 34) example -------------------------------------------------

test_that("LTabr works on UN 1982 (p. 34) example", {
    
    Mx <- c(.23669, .04672, .00982, .00511, .00697, .01036, .01169,
            .01332, .01528, .01757, .02092, .02517, .03225, .04241, 
            .06056, .08574, .11840, .16226, .23745)
    
    excheckUN <- c(35.000, 42.901, 47.190, 44.438, 40.523, 36.868, 
                   33.691, 30.567, 27.500, 24.485, 21.504, 18.599,
                   15.758, 13.080, 10.584, 8.466, 6.729, 5.312, 4.211)
    
    AgeInt <- inferAgeIntAbr(vec = Mx)
    
    # generate two variants: with and without PG's variants for ages 5-14
    UNLT1 <- LTabr(nMx = Mx,
                   Age = c(0,1,seq(5,85,by=5)),
                   AgeInt = AgeInt,
                   axmethod = "un",
                   Sex = "m",
                   mod = FALSE)
    
    UNLT2 <- LTabr(nMx = Mx,
                   Age = c(0,1,seq(5,85,by=5)),
                   AgeInt = AgeInt,
                   axmethod = "un",
                   Sex = "m",
                   mod = TRUE)
    
    # test
    expect_equal(
        UNLT1$ex,
        excheckUN,
        tolerance = .01
    )
})



# Mrotpack example --------------------------------------------------------

test_that("LTabr works on Mortpak exapmle (United Nations 1988, p. 82)", {
    
    MPnMx <- c(0.12846, 0.02477, 0.00603, 0.0034,
               0.00417, 0.00513, 0.00581, 0.00645, 0.00725,
               0.00813, 0.00913, 0.01199, 0.01647,
               0.0256, 0.04047, 0.06624, 0.10638, 0.19611)
    
    MPexcheck <- c(49.997, 55.675, 57.245, 53.921,
                   49.803, 45.799, 41.922, 38.084, 34.249,
                   30.420, 26.578, 22.701, 18.945,
                   15.349, 12.095, 9.240, 6.903, 5.099)
    
    # First with lifetable extention to 100
    MP_UNLT100 <- LTabr(
        nMx = MPnMx,
        Age = c(0, 1, seq(5, 80, by = 5)),
        AgeInt = inferAgeIntAbr(vec = MPnMx),
        axmethod = "un",
        Sex = "f",
        mod = FALSE,
        OAnew = 100
    )
    
    # lifetable to original open age group
    MP_UNLT80 <- LTabr(
        nMx = MPnMx,
        Age = c(0, 1, seq(5, 80, by = 5)),
        AgeInt = inferAgeIntAbr(vec = MPnMx),
        axmethod = "un",
        Sex = "f",
        mod = FALSE,
        OAnew = 80
    )
    
    # same, but truncated at 60
    MP_UNLT60 <- LTabr(
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
    expect_equal(
        MP_UNLT100$ex[1:length(MPexcheck)],
        MPexcheck,
        tolerance = .002
    )
    
    expect_equal(
        MP_UNLT80$ex,
        MPexcheck,
        tolerance = .002
    )
    
    # identical results irrespective of max age
    expect_equal(
        MP_UNLT60$ex,
        MP_UNLT80$ex[1:length(MP_UNLT60$ex)],
        tolerance = 1e-12
    )
    
    expect_equal(
        MP_UNLT80$ex,
        MP_UNLT100$ex[1:length(MP_UNLT80$ex)],
        tolerance = 1e-12
    )
})
