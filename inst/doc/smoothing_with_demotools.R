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
     ylab = 'Porportion of population',
     xlab = 'Single age',
     main = 'Unsmoothed population') 

  

## ----out.width='\\textwidth', fig.width= 6, fig.height=6, fig.align='center'----
cf <- agesmth(Value = Value, 
	Age = Age, 
	method = "Carrier-Farrag", 
	OAG = TRUE)

plot(seq(0,100,5),cf/sum(cf), type= 'l',
     ylab = 'Porportion of population',
     xlab = 'Age-group',main = 'Smoothed population with Carrier-Farrag',xaxt='n')
axis(1, labels = paste0(seq(0,100,5),'-',seq(4,104,5)), at =seq(0,100,5))
  

## ----out.width='\\textwidth', fig.width= 6, fig.height=6, fig.align='center'----
cf <- agesmth(Value = Value, 
	Age = Age, 
	method = "Arriaga", 
	OAG = TRUE)

plot(seq(0,100,5),cf/sum(cf), type= 'l',
     ylab = 'Porportion of population',
     xlab = 'Age-group',main = 'Smoothed population with Arriaga',xaxt='n')
axis(1, labels = paste0(seq(0,100,5),'-',seq(4,104,5)), at =seq(0,100,5))
  

## ----out.width='\\textwidth', fig.width= 6, fig.height=6, fig.align='center'----
cf <- agesmth(Value = Value, 
	Age = Age, 
	method = "Karup-King-Newton", 
	OAG = TRUE)

plot(seq(0,100,5),cf/sum(cf), type= 'l',
     ylab = 'Porportion of population',
     xlab = 'Age-group',main = 'Smoothed population with Karup-King-Newton ',xaxt='n')
axis(1, labels = paste0(seq(0,100,5),'-',seq(4,104,5)), at =seq(0,100,5))
  

