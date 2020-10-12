
# Author: ilya
###############################################################################
context("test-ageheap")

# PASEX spreadsheet SINGAGE

test_that("check_heaping_whipple replicates SINGAGE males", {
    
    Age <- 0:99

    # test
    expect_equal(
        check_heaping_whipple(
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
test_that("is_single() applied only to tested age range",{

# should induce perfect equality also
df <- data.frame(
		age = c(0:90, 100),
		count = 10
)

expect_equal(check_heaping_whipple(df$count, df$age, ageMin = 25, ageMax = 65),1)
expect_equal(check_heaping_myers(df$count, df$age, ageMin = 25, ageMax = 65),0)
expect_equal(check_heaping_bachi(df$count, df$age, ageMin = 25, ageMax = 65),0)
expect_equal(check_heaping_coale_li(df$count, df$age, ageMin = 60, ageMax = 89),1)
expect_equal(check_heaping_noumbissi(df$count, df$age, ageMin = 25, ageMax = 65),1)
expect_equal(check_heaping_spoorenberg(df$count, df$age, ageMin = 25, ageMax = 65),0)

# should induce errors
# no longer throw errors because these are lowered where necessary
# expect_error(check_heaping_whipple(df$count, df$age, ageMin = 25, ageMax = 100))
# expect_error(check_heaping_myers(df$count, df$age, ageMin = 25, ageMax = 100))
# expect_error(check_heaping_bachi(df$count, df$age, ageMin = 25, ageMax = 100))
# expect_error(check_heaping_coale_li(df$count, df$age, ageMin = 60, ageMax = 100))
# expect_error(check_heaping_noumbissi(df$count, df$age, ageMin = 25, ageMax = 100))
# expect_error(check_heaping_spoorenberg(df$count, df$age, ageMin = 25, ageMax = 100))
})