
# Author: tim
###############################################################################

library(devtools)
document("/home/tim/git/DemoTools")

load_all("/home/tim/git/DemoTools")
use_appveyor("/home/tim/git/DemoTools")
use_travis("/home/tim/git/DemoTools")

use_coverage(pkg = "/home/tim/git/DemoTools", type = c("codecov", "coveralls"))

check("/home/tim/git/DemoTools")
