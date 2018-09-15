# Author: ilya
###############################################################################
context("test-zigzag")

# Feeney, G. 2013 "Removing "Zigzag" from Age Data," http://demographer.com/white-papers/2013-removing-zigzag-from-age-data/

test_that("zigzag works",{
    
    Value <- c(13331, 4151, 1746, 1585, 3859, 8354, 11146, 12076,
               12216, 12016, 12473, 11513, 12899, 11413, 12710, 
               11516, 11408, 6733, 4031, 2069)
    
    Age <- c(0, 1, seq(5, 90, by = 5))

    # check result using results frozen in Feeney spreadsheet
    # after fixing probable cell range 'error'
    p.feeney <- c(0.0235802695087692, 0.0724286618207911,
                  0.0242327829742267, 0.0883411499065237)
    
    ans      <- 106.1147411629

    
    expect_equal(
        zigzag_min(Value, Age, 40,80,p.feeney),
        ans,
        tolerance = 1e-6
    )
})

