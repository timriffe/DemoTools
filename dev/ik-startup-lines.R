# Author: IK
# Date: 2018-10-12
# Sturtup lines to work on package development
################################################################################

library(magrittr)
library(devtools)
library(testthat)
# rebuild the latest version of DemoTools
install_github("timriffe/DemoTools")
library(DemoTools)