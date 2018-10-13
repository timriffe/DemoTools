# Author: ilya
###############################################################################
context("test-grabill")

test_that("grabillExpand works",{
    expect_equal(grabillExpand(pop5_mat, OAG = TRUE)[1:3, 1], 
                 c(0.0826, 0.0903, 0.0932), tolerance = 1e-3)
    expect_equal(grabillExpand(pop5_mat, OAG = TRUE)[101, 21], 1)
})


test_that("grabill works",{
    popmat <- matrix(pop1m_ind, dimnames = list(0:100, NULL))
    expect_equal(as.vector(head(grabill(popmat), 2)), 
                 c(9350426.61, 10050060.85), 
                 tolerance = 1e-2)
})


