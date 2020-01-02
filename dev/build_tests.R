# testing dependencies and so forth

# -----------------------------------------
# visualize individual function dependencies in DemoTools
#install.packages("yaml")
#install.packages("visNetwork")
#devtools::install_github("datastorm-open/DependenciesGraphs")
library(DependenciesGraphs)

library(DemoTools) # The package we want to explore
# before tinkering with an older function, note which functions depend on it
deps <- funDependencies("package:DemoTools","inferAgeIntAbr")
plot(deps)

# What about our own package dependency tree?
deps <- tools::package_dependencies(
  packages="DemoTools",
  which = c("Depends","Imports"),
  recursive = TRUE, 
  db = installed.packages()) 
length(deps[[1]])  # 112 !!

# how many deps causes by each first order dep of DemoTools?
d1 <- tools::package_dependencies(
  packages="DemoTools",
  which = c("Depends","Imports"),
  recursive = FALSE, 
  db = installed.packages()) 
library(magrittr)

d2 <- tools::package_dependencies(
  packages=d1$DemoTools,
  which = c("Depends","Imports"),
  recursive = TRUE, 
  db = installed.packages())

dcoutns <- d2 %>% 
  lapply(length) %>%   
  unlist() %>% 
  sort() 

barplot(dcounts, horiz = TRUE, las = 1)

d2$ungroup

# Counterfactuals: eliminiting which package would have the biggest payoff?
n <- length(d2)
N <- length(deps[[1]])
eliminable <- rep(0,n)
for (i in 1:n){
  eliminable[i] <-
    unlist(d2[-i]) %>% unique() %>% length()
}
str(d2)

names(eliminable) <- names(d2)
N - eliminable
unlist(d2) %>% unique()