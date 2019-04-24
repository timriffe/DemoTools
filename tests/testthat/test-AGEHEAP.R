
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
test_that("is_single() applied only to tested age range",{

# should induce perfect equality also
df <- data.frame(
		age = c(0:90, 100),
		count = 10
)

expect_equal(Whipple(df$count, df$age, ageMin = 25, ageMax = 65),1)
expect_equal(Myers(df$count, df$age, ageMin = 25, ageMax = 65),0)
expect_equal(Bachi(df$count, df$age, ageMin = 25, ageMax = 65),0)
expect_equal(CoaleLi(df$count, df$age, ageMin = 60, ageMax = 89),1)
expect_equal(Noumbissi(df$count, df$age, ageMin = 25, ageMax = 65),1)
expect_equal(Spoorenberg(df$count, df$age, ageMin = 25, ageMax = 65),0)

# should indice errors
expect_error(Whipple(df$count, df$age, ageMin = 25, ageMax = 100))
expect_error(Myers(df$count, df$age, ageMin = 25, ageMax = 100))
expect_error(Bachi(df$count, df$age, ageMin = 25, ageMax = 100))
expect_error(CoaleLi(df$count, df$age, ageMin = 60, ageMax = 100))
expect_error(Noumbissi(df$count, df$age, ageMin = 25, ageMax = 100))
expect_error(Spoorenberg(df$count, df$age, ageMin = 25, ageMax = 100))
})