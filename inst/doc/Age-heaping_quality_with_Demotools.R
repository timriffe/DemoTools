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

## ------------------------------------------------------------------------
P5 <- groupAges(pop1m_ind, N = 5, 0:100)
A5 <- seq(0, 100, by = 5)
plot(A5,P5,type='o')

## ------------------------------------------------------------------------
zero_pref_sawtooth(pop1m_ind, Age = 0:100, ageMin = 40, ageMax = 80)

## ------------------------------------------------------------------------
r5 <- five_year_roughness(pop1m_ind, Age = 0:100, ageMin = 40, ageMax = 80)
r5

