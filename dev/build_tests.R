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

# what file is each function in?

A <- lsf.str("package:DemoTools")
devtools::load_all()
attr(attr(inferAgeIntAbr,"srcref"),"srcfile")


get_file <- function(filename_no_quotes){
  attr(attr(filename_no_quotes,"srcref"),"srcfile")
}
get_file(inferAgeIntAbr)
lapply(ls("package:DemoTools"),function(x) eval() %>% is.function)
is.function(DemoTools::`:=`)

# list of functions by script

scripts  <- dir("R")
files <- list()

for (i in scripts){
test.env <- new.env()
sys.source(paste0("R/",i), envir = test.env)
A <- lsf.str(envir=test.env)
funs <- lapply(A,"[[",1) %>% unlist()
files[[i]] <- funs
rm(A)
rm(test.env)
}
library(tidyverse)

lengths <- lapply(files, length) %>% unlist()
DF      <- tibble(Script = rep(names(lengths), times = lengths),
                  Function = unlist(files))
DF %>% 
  arrange(Function) %>% 
  write_csv("dev/FunctionInventory.csv")
