## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----out.width='\\textwidth', fig.width= 6, fig.height=6, fig.align='center'----
library(DemoTools)

# 'pop1m_ind' available as package data
Value <- pop1m_ind
Age   <- 0:(length(Value)-1)
  
plot(Age, Value/sum(Value), type = 'l',
     ylab = 'Proportion of population',
     xlab = 'Single age',
     main = 'Unsmoothed population') 

  

## ----out.width='\\textwidth', fig.width= 6, fig.height=6, fig.align='center'----
cf <- agesmth(Value = Value, 
	Age = Age, 
	method = "Carrier-Farrag", 
	OAG = TRUE)

plot(seq(0,100,5),cf/sum(cf), type= 'l',
     ylab = 'Proportion of population',
     xlab = 'Age-group',main = 'Smoothed population with Carrier-Farrag',xaxt='n')
axis(1, labels = paste0(seq(0,100,5),'-',seq(4,104,5)), at =seq(0,100,5))
  

## ----out.width='\\textwidth', fig.width= 6, fig.height=6, fig.align='center'----
cf <- agesmth(Value = Value, 
	Age = Age, 
	method = "Arriaga", 
	OAG = TRUE)

plot(seq(0,100,5),cf/sum(cf), type= 'l',
     ylab = 'Proportion of population',
     xlab = 'Age-group',main = 'Smoothed population with Arriaga',xaxt='n')
axis(1, labels = paste0(seq(0,100,5),'-',seq(4,104,5)), at =seq(0,100,5))
  

## ----out.width='\\textwidth', fig.width= 6, fig.height=6, fig.align='center'----
cf <- agesmth(Value = Value, 
	Age = Age, 
	method = "Karup-King-Newton", 
	OAG = TRUE)

plot(seq(0,100,5),cf/sum(cf), type= 'l',
     ylab = 'Proportion of population',
     xlab = 'Age-group',main = 'Smoothed population with Karup-King-Newton ',xaxt='n')
axis(1, labels = paste0(seq(0,100,5),'-',seq(4,104,5)), at =seq(0,100,5))
  

## ----out.width='\\textwidth', fig.width= 6, fig.height=6, fig.align='center'----
 un_result <- united_nations_smth(Value = Value,Age = Age,OAG = T)
plot(seq(0,100,5),un_result/sum(un_result,na.rm = T), type= 'l',
    ylab = 'Proportion of population',
     xlab = 'Age-group',main = 'Smoothed population with UN Method',xaxt='n')
axis(1, labels = paste0(seq(0,100,5),'-',seq(4,104,5)), at =seq(0,100,5))

## ----out.width='\\textwidth', fig.width= 6, fig.height=6, fig.align='center'----
 strong_result <- strong_smth(Value = Value,Age = Age,OAG = T)
plot(seq(0,100,5),strong_result/sum(strong_result,na.rm = T), type= 'l',
    ylab = 'Proportion of population',
     xlab = 'Age-group',main = 'Smoothed population with Strong formula',xaxt='n')
axis(1, labels = paste0(seq(0,100,5),'-',seq(4,104,5)), at =seq(0,100,5))

## ----out.width='\\textwidth', fig.width= 6, fig.height=6, fig.align='center'----
 zz_result <- zigzag_smth(Value, Age, OAG = TRUE, ageMin = 0, ageMax = 100)
plot(seq(0,100,5),strong_result/sum(strong_result,na.rm = T), type= 'l',
    ylab = 'Proportion of population',
     xlab = 'Age-group',main = 'Smoothed population with Zig Zag formula',xaxt='n')
axis(1, labels = paste0(seq(0,100,5),'-',seq(4,104,5)), at =seq(0,100,5))

## ----out.width='\\textwidth', fig.width= 6, fig.height=6, fig.align='center'----
 poly_result <- agesmth1(Value, Age, method="poly", OAG = T)
plot(0:100,poly_result, type= 'l',
    ylab = 'Proportion of population',
     xlab = 'Age-group',main = 'Smoothed population with Poly formula',xaxt='n')
points(0:100,Value)
axis(1, labels =seq(0,100,5), at =seq(0,100,5))

## ----out.width='\\textwidth', fig.width= 6, fig.height=6, fig.align='center'----
 loess_result <- agesmth1(Value, Age, method="loess", OAG = T)
plot(0:100,loess_result, type= 'l',
    ylab = 'Proportion of population',
     xlab = 'Age-group',main = 'Smoothed population with LOESS formula',xaxt='n')
points(0:100,Value)
axis(1, labels =seq(0,100,5), at =seq(0,100,5))

