
# Author: tim
###############################################################################

library(devtools)
# do this whenever new functions are added to /R, or whenever roxygen is updated
document("/home/tim/git/DemoTools")

# run this to get access to already-written functions
load_all("/home/tim/git/DemoTools")

# do this whenever major changes happen
check("/home/tim/git/DemoTools")


# these created the necessary files to run automatic remote code testing
#use_appveyor("/home/tim/git/DemoTools")
#use_travis("/home/tim/git/DemoTools")
#use_coverage(pkg = "/home/tim/git/DemoTools", type = c("codecov", "coveralls"))
