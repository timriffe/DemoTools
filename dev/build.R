
# Author: tim
###############################################################################
# old:
#  rstan (>= 2.18.1),
# new:
#  rstan (>= 2.26.1),
# Tried adding this to Remotes:
# list in DESCRIPTION:
# github::hsbadr/rstan/StanHeaders@develop,
# github::hsbadr/rstan/rstan/rstan@develop

shhh <- function(expr){
	capture.output(x <- suppressPackageStartupMessages(
					suppressMessages(suppressWarnings(expr))))
	invisible(x)
}

library(devtools)
library(TimUtils)
library(magrittr)
library(git2r)
#install.packages("backports")
#install.packages("roxygen2")
#install_github("hadley/devtools")
#install_github("timriffe/TimUtils")

# do this whenever new functions are added to /R, or whenever roxygen is updated
devtools::document()

# do this whenever the vignette text is updated
devtools::build_vignettes()

# devtools::install_github("r-lib/pkgdown")
usethis::proj_activate(here::here())
pkgdown::build_site()
  
versionIncrement(
		major = FALSE,       # only for releases
		mid = FALSE,         # major functionality added
		minor = TRUE,        # whenever documentation renewed, any patch, tweak, or fix
		maxdigits = c(2,2,3),# maybe 4 required?
		README = TRUE)       # update README dev version badge

# add line to immediately commit and tag.

# D <- readLines("DESCRIPTION") 
# vs <- D[grepl(D,pattern = "Version: ")]  %>%  gsub(pattern = "Version: ", replacement = "")  %>% 
#   paste0("v",.)
# commit(message = vs,all=TRUE)
# tag(name =vs,message = vs)
# push(refspec = vs)

# https://raw.githubusercontent.com/timriffe/DemoTools/59a0f4e50b7696c185a3c9d4e582426f88aac84f/DESCRIPTION


# run this to get access to already-written functions
shhh(load_all())

# usethis::use_build_ignore(c("docs"))
# do this whenever major changes happen
Sys.setenv('_R_CHECK_SYSTEM_CLOCK_' = 0)
devtools::check(force_suggests = TRUE)


source("version_lookup.R")
update_lookup()
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


