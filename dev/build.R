
# Author: tim
###############################################################################
me <- system("whoami",TRUE)
if (me == "tim"){
	setwd("/home/tim/git/DemoTools")
}
# add in your own statement like this if you like


shhh <- function(expr){
	capture.output(x <- suppressPackageStartupMessages(
					suppressMessages(suppressWarnings(expr))))
	invisible(x)
}

library(devtools)
library(TimUtils)
#install_github("hadley/devtools")
#install_github("timriffe/TimUtils/TimUtils")
# do this whenever new functions are added to /R, or whenever roxygen is updated
devtools::document()
versionIncrement(
		major = FALSE,       # only for releases
		mid = FALSE,         # major functionality added
		minor = TRUE,        # whenever documentation renewed, any patch, tweak, or fix
		maxdigits = c(2,2,3),# maybe 4 required?
		README = TRUE)       # update README dev version badge

# run this to get access to already-written functions
shhh(load_all())

# do this whenever major changes happen
devtools::check()
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
length(dir("/home/tim/git/DemoTools/man"))

#library(badger)
#badge_devel("timriffe/DemoTools", "yellow")

versionIncrement(
		major = FALSE,       # only for releases
		mid = TRUE,         # major functionality added
		minor = FALSE,        # whenever documentation renewed, any patch, tweak, or fix
		maxdigits = c(2,2,3),# maybe 4 required?
		README = TRUE)  