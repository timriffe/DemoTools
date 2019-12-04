# Author: ilya
###############################################################################
context("test-grabill")

test_that("grabillExpand works",{
  a5 <- as.integer(rownames(pop5_mat))
    expect_equal(graduate_grabill_expand(pop5_mat[,1], Age = a5, OAG = TRUE)[1:3, 1], 
                 c(0.0826, 0.0903, 0.0932), tolerance = 1e-3)
    expect_equal(graduate_grabill_expand(pop5_mat[,1], Age = a5, OAG = TRUE)[101, 21], 1)
})


test_that("grabill works",{
    check <- c(9350425.61, 10050059.85)
    names(check) <- c(0,1)
    expect_equal(c(graduate_grabill(pop1m_ind, 0:100, TRUE)[1:2]), 
                 check, 
                 tolerance = 1e-2)
})
