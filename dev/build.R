
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
#install_github("hadley/devtools")
# do this whenever new functions are added to /R, or whenever roxygen is updated
devtools::document()

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