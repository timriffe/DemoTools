
# Author: tim
###############################################################################

library(devtools)
document("/home/tim/git/DemoTools/DemoTools")

load_all("/home/tim/git/DemoTools/DemoTools")
use_appveyor("/home/tim/git/DemoTools/DemoTools")
use_travis("/home/tim/git/DemoTools/DemoTools")

use_coverage(pkg = "/home/tim/git/DemoTools/DemoTools", type = c("codecov", "coveralls"))

check("/home/tim/git/DemoTools/DemoTools")
?devtools::use_appveyor