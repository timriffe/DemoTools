
# Author: tim
###############################################################################

library(devtools)
#install_github("hadley/devtools")
# do this whenever new functions are added to /R, or whenever roxygen is updated
document("/home/tim/git/DemoTools")

# run this to get access to already-written functions
load_all("/home/tim/git/DemoTools")

# do this whenever major changes happen
check("/home/tim/git/DemoTools")
#build(pkg = "/home/tim/git/DemoTools", path = "/home/tim/Desktop")
#?devtools::build
#devtools::use_testthat("/home/tim/git/DemoTools")
#dir("/home/tim/git/DemoTools/man")

# these created the necessary files to run automatic remote code testing
#use_appveyor("/home/tim/git/DemoTools")
#use_travis("/home/tim/git/DemoTools")
#use_coverage(pkg = "/home/tim/git/DemoTools", type = c("codecov", "coveralls"))



# how many functions are there?
# length(dir("/home/tim/git/DemoTools/man"))