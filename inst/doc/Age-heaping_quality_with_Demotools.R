## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6, 
  fig.height = 6,
  fig.align = "center"
)


## ------------------------------------------------------------------------
library(DemoTools)

# 'pop1m_ind' available as package data
Age   <- 0:(length(pop1m_ind)-1)
 
plot(Age, pop1m_ind, type = 'o')  
  

## ------------------------------------------------------------------------
Wi <-  Whipple(pop1m_ind, Age, ageMin = 25, ageMax = 60, digit = c(0, 5))
Wi

## ------------------------------------------------------------------------
Mi <- Myers(pop1m_ind, Age, ageMin = 20, ageMax = 90) 
Mi


## ------------------------------------------------------------------------
Bi <- Bachi(pop1m_ind, Age, ageMin = 20, ageMax = 90)
Bi

## ------------------------------------------------------------------------
Ni <- Noumbissi(pop1m_ind, Age, digit = 0) 
Ni

## ------------------------------------------------------------------------
Si <- Spoorenberg(pop1m_ind, Age, ageMin = 20, ageMax = 64)
Si

## ------------------------------------------------------------------------
CLi <- CoaleLi(pop1m_ind, Age, ageMin = 60, ageMax = max(Age), terms = 5, digit = 0)
CLi

## ------------------------------------------------------------------------
Ji <- Jdanov(pop1m_ind, Age, Agei = c(95,100,105))
Ji

