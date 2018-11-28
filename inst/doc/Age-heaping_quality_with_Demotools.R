## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  tidy = FALSE,
  comment = "#>",
  fig.width = 6, fig.height = 6
)

## ------------------------------------------------------------------------
library(DemoTools)

# 'pop1m_ind' available as package data
Value <- pop1m_ind
Age   <- 0:(length(Value)-1)
  
plot(Age, Value, type = 'l')  
  

## ------------------------------------------------------------------------
Whipple <-  Whipple(Value, Age, ageMin = 23, ageMax = 62, digit = c(0,5))
Whipple

## ------------------------------------------------------------------------
Myers <- Myers(Value, Age, ageMin = 10, ageMax = 90) 
Myers


## ------------------------------------------------------------------------
Bachi <- Bachi(Value, Age, ageMin = 10, ageMax = 90)

## ------------------------------------------------------------------------
CoaleLi <- CoaleLi(Value, Age, ageMin = 60, ageMax = max(Age), terms = 5, digit = 0)
CoaleLi

## ------------------------------------------------------------------------
Noumbissi <- Noumbissi(Value, Age, digit = 0) 
Noumbissi

## ------------------------------------------------------------------------
Spoorenberg <- Spoorenberg(Value, Age, ageMin = 20, ageMax = 64)
Spoorenberg

