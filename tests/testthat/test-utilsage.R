# Author: ilya
###############################################################################
context("test-utilsage")



test_that("age2ageN works",{
    age1 <- seq(0,100,by=5)
    (ageN1 <- age2ageN(age1, OAG = FALSE))
    (ageN2 <- age2ageN(age1, OAG = TRUE))
    
    expect_equal(length(ageN1),  105)
    expect_equal(length(ageN2),  101)
})


test_that("int2ageN works",{
    int5 <- rep(5,21)
    (ageN1 <- int2ageN(int5, OAG = FALSE))
    (ageN2 <- int2ageN(int5, OAG = TRUE))
    
    expect_equal(length(ageN1),  105)
    expect_equal(length(ageN2),  101)
})



test_that("calcAgeAbr works",{
    Age <- 0:70
    ageN1 <- calcAgeAbr(Age)

    expect_equal(length(ageN1),  71)
    expect_equivalent(
        as.vector(head(table(ageN1), 3)),
        c(1, 4, 5)
    )
})


test_that("calcAgeN works",{
    Age <- 0:70
    ageN1 <- calcAgeN(Age,5,0)
    
    expect_equal(length(ageN1),  71)
    expect_equivalent(
        as.vector(head(table(ageN1), 3)),
        c(5, 5, 5)
    )
})


test_that("inferAgeIntAbr works",{
    vec <- runif(20)
    one <- inferAgeIntAbr(vec = vec)
    two <- inferAgeIntAbr(vec = vec, OAG = TRUE)
    
    expect_equivalent(head(one, 3), c(1, 4, 5))
    expect_true(is.na(tail(two, 1)))

})


test_that("maxA2abridged works",{
    expect_true(all(maxA2abridged(100) == maxA2abridged(102)))
})


test_that("int2age works",{
    AgeInt <- c(1, 4, rep(5, 17))
    one <- int2age(AgeInt)
    expect_gt(tail(one, 1), head(one, 1))
    expect_length(one, 19)
})



test_that("groupAges works",{
    Age <- 0:100
    expect_equivalent(
        groupAges(pop1m_ind, N = 5) [1:3],
        c(52360652, 57419164, 51947630)
    )
    expect_equivalent(
        groupAges(pop1m_ind, N = 5, shiftdown = 1) [1:3],
        c(40488149, 61180144, 51043772)
    )
    expect_equivalent(
        groupAges(pop1m_ind, N = 5, shiftdown = 2) [1:3],
        c(28606305, 58749766, 57901419)
    )
    expect_equivalent(
        groupAges(pop1m_ind, N = 5, shiftdown = 3) [1:3],
        c(17016196, 60305957, 54630220)
    )
    expect_equivalent(
        groupAges(pop1m_ind, N = 5, shiftdown = 4) [1:3],
        c(9544406, 55784596, 59761861)
    )
    expect_equivalent(
        tail(groupAges(pop1m_ind, N = 5, OAnew = 80), 3),
        c(5535950, 2102284, 3324624)
    )
})



test_that("is_single works",{
    Age <- 0:99
    Age2 <- c(0:10, 0:10)
    Age3 <- seq(0, 80, by = 5)
    Age4 <- seq(0, 10, by = .5)
    
    expect_true(is_single(Age))
    expect_false(is_single(Age2))
    expect_false(is_single(Age3))
    expect_false(is_single(Age4))
})


test_that("is_abridged works",{
    vec <- c(0, 1, 5, 10, 15, 20, 25)
    expect_true(is_abridged(vec))
    expect_false(is_abridged(vec[-2]))
    expect_true(is_abridged(tail(vec, -1)))
    expect_true(is_abridged(tail(vec, -2)))
    expect_false(is_abridged(vec[-5]))
})



test_that("rescaleAgeGroups works",{
    # just to make a point about arbitrary integer age widths in both pop1 and pop2
    # note if pop1 is in single ages and pop2 is in groups things work much cleaner.
    set.seed(3)
    AgeIntRandom <- sample(1:5, size = 15,replace = TRUE)
    AgeInt5      <- rep(5, 9)
    original     <- runif(45, min = 0, max = 100)
    pop1         <- groupAges(original, 0:45, AgeN = int2ageN(AgeIntRandom, FALSE))
    pop2         <- groupAges(original, 0:45, AgeN = int2ageN(AgeInt5, FALSE))
    # inflate (in this case) pop2
    perturb      <- runif(length(pop2), min = 1.05, max = 1.2)
    
    pop2         <- pop2 * perturb
    
    # a recursively constrained solution
    pop1resc <- rescaleAgeGroups(Value1 = pop1,
                                  AgeInt1 = AgeIntRandom,
                                  Value2 = pop2,
                                  AgeInt2 = AgeInt5,
                                  splitfun = splitUniform,
                                  recursive = TRUE)
    # a single pass adjustment (no recursion)
    pop1resc1 <- rescaleAgeGroups(Value1 = pop1,
                                   AgeInt1 = AgeIntRandom,
                                   Value2 = pop2,
                                   AgeInt2 = AgeInt5,
                                   splitfun = splitUniform,
                                   recursive = FALSE)
    
    newN        <- splitUniform(pop1resc, AgeInt = AgeIntRandom)
    AgeS        <- names2age(newN)
    AgeN2       <- rep(int2age(AgeInt5), times = AgeInt5)
    check       <- groupAges(newN, AgeS, AgeN = AgeN2)
    
    expect_equal(check, pop2, tolerance = 1e-3)
    expect_equal(sum(pop2), sum(pop1resc1))
    expect_equal(sum(pop1resc), sum(pop1resc1))
})


