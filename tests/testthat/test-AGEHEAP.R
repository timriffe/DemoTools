
# Author: ilya
###############################################################################
context("test-ageheap")

# PASEX spreadsheet SINGAGE

test_that("Whipple replicates SINGAGE males", {
    
    Age <- 0:99

    # test
    expect_equal(
        Whipple(
            Value = pop1m_pasex, 
            Age = Age, 
            ageMin = 25, 
            ageMax = 60, 
            digit = c(0, 5)
        ), 
        2.34, 
        tolerance = .01
    )

})
