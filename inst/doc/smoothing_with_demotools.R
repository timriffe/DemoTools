## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(DemoTools)

# 'pop1m_ind' available as package data
Value <- pop1m_ind
Age   <- 0:(length(Value)-1)
  
plot(Age, Value, type = 'l')  
  

