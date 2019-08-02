# Author: ilya

context("test-agesmth")


# carrier_farrag_smth -----------------------------------------------------

test_that("carrier_farrag_smth works",{
    
    Ages         <- seq(0, 80, by = 5)
    
    CFmales <- carrier_farrag_smth(pop5m_pasex, Ages, OAG = TRUE)
    
    CFtest <- c(NA, NA, 346290, 287083, 285855, 261082, 237937,
                202809, 162973, 125720, 88730, 67352, 55187, 40657, NA, NA, NA)
    
    # test
    expect_equivalent(CFmales, CFtest, tolerance = 1)
    expect_equal(CFmales %>% sum, CFtest %>% sum, tolerance = 1)
    expect_true(all(CFmales > 0, na.rm = T))

})



# kkn_smth ----------------------------------------------------------------

test_that("kkn_smth works",{
    
    Ages <- seq(0, 80, by = 5)
    
    KKNtest <- c(NA, NA, 354871, 278502, 285508, 261429, 236513 ,
                 204233, 162138, 126555, 90094, 65988, 54803, 41041, NA, NA, NA)
    
    KKNmales <- kkn_smth(pop5m_pasex, Ages, OAG = TRUE)
    
    # test
    expect_equivalent(KKNmales, KKNtest, tolerance = 1)
    expect_equal(KKNmales %>% sum, KKNtest %>% sum, tolerance = 1)
    expect_true(all(KKNmales > 0, na.rm = T))
    
})


# arriaga_smth ----------------------------------------------------------------

test_that("arriaga_smth works",{
    
    Ages <- seq(0, 80, by = 5)
    
    Amales <- arriaga_smth(Value = pop5m_pasex, Age = Ages, OAG = TRUE)
    
    #' # PAS spreadsheet result:
    Atest <- c(662761, 495126, 345744, 287629, 285919, 261018, 237469, 203277,
    161733, 126960, 88586, 67496, 54587, 41257, 28790, 17189, 34729)
    
    # test
    expect_equivalent(Amales, Atest, tolerance = 1)
    expect_equal(Amales %>% sum, Atest %>% sum, tolerance = 1)
    expect_true(all(Amales > 0, na.rm = T))
    
})


# strong_smth ----------------------------------------------------------------

test_that("strong_smth works",{
    
    Ages <- seq(0, 80, by = 5)
    
    Stest <- c(646617, 511270, 386889, 317345, 273736, 240058, 218645, 188297, 
               153931, 124347, 93254, 71858, 53594, 39721, 27887, 18092, 34729)
    
    Smales <- strong_smth(pop5m_pasex, Ages, OAG = TRUE)
    
    # test
    expect_equivalent(Smales, Stest, tolerance = 1)
    expect_equal(Smales %>% sum, Stest %>% sum, tolerance = 1)
    expect_true(all(Smales > 0, na.rm = T))
    
})


#' Age <- c(0,1,seq(5,90,by=5))
#' # defaults
#' zz <- zigzag_smth(dth5_zigzag, Age, OAG = TRUE, ageMin = 40, ageMax = 90)