
# Author: tim
###############################################################################

shhh <- function(expr){
	capture.output(x <- suppressPackageStartupMessages(
					suppressMessages(suppressWarnings(expr))))
	invisible(x)
}

library(devtools)
library(TimUtils)
#install.packages("backports")
#install.packages("roxygen2")
#install_github("hadley/devtools")
#install_github("timriffe/TimUtils/TimUtils")

# do this whenever new functions are added to /R, or whenever roxygen is updated
devtools::document()
# do this whenever the vignette text is updated
devtools::build_vignettes()

# devtools::install_github("r-lib/pkgdown")
pkgdown::build_site()
  
versionIncrement(
		major = FALSE,       # only for releases
		mid = FALSE,         # major functionality added
		minor = TRUE,        # whenever documentation renewed, any patch, tweak, or fix
		maxdigits = c(2,2,3),# maybe 4 required?
		README = TRUE)       # update README dev version badge

# run this to get access to already-written functions
shhh(load_all())

# usethis::use_build_ignore(c("docs"))
# do this whenever major changes happen
devtools::check(force_suggests = TRUE, manual = FALSE)

  #build(pkg = "/home/tim/git/DemoTools", path = "/home/tim/Desktop")
#?devtools::build
#devtools::use_testthat("/home/tim/git/DemoTools")
#dir("/home/tim/git/DemoTools/man")
install_github("timriffe/DemoTools")
# these created the necessary files to run automatic remote code testing
#use_appveyor("/home/tim/git/DemoTools")
#use_travis("/home/tim/git/DemoTools")
#use_coverage(pkg = "/home/tim/git/DemoTools", type = c("codecov", "coveralls"))

# to update statement in README.md: approximates nr functions avail
length(dir(here::here("man")))

#library(badger)
#badge_devel("timriffe/DemoTools", "yellow")


# mid-level increment (if a new top-level function is added)
versionIncrement(
		major = FALSE,        # only for releases
		mid = TRUE,           # major functionality added
		minor = FALSE,        # whenever documentation renewed, any patch, tweak, or fix
		maxdigits = c(2,2,3), # maybe 4 required?
		README = TRUE)  

# top level increment
versionIncrement(
		major = TRUE,         # only for releases
		mid = FALSE,          # major functionality added
		minor = FALSE,        # whenever documentation renewed, any patch, tweak, or fix
		maxdigits = c(2,2,3), # maybe 4 required?
		README = TRUE)  
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
i<- 1
# for setting options
#candidates <- c( Sys.getenv("R_PROFILE"),
#		file.path(Sys.getenv("R_HOME"), "etc", "Rprofile.site"),
#		Sys.getenv("R_PROFILE_USER"),
#		file.path(getwd(), ".Rprofile") )
#
#Filter(file.exists, candidates)

# NOTE TO SELF
# try goodpractice package

# Extra once-off checks

# checks run Aug 13, 2018
check_win_devel()      # OK
check_win_release()    # OK
check_win_oldrelease() # OK

check_rhub(email = "tim.riffe@gmail.com", interactive = FALSE)  # sent

library(spelling)
spell_check()

# for notable moments of stability and cleanliness:
release()


