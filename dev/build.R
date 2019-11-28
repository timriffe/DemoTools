
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


# do this whenever major changes happen
devtools::check(force_suggests=TRUE,manual=FALSE)

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
# visualize function dependencies in DemoTools
#install.packages("yaml")
#install.packages("visNetwork")
#devtools::install_github("datastorm-open/DependenciesGraphs")
library(DependenciesGraphs)

library(DemoTools) # The package we want to explore
# before tinkering with an older function, note which functions depend on it
deps <- funDependencies("package:DemoTools","mxax2qx_Backstop")
plot(deps)


# for setting options
#candidates <- c( Sys.getenv("R_PROFILE"),
#		file.path(Sys.getenv("R_HOME"), "etc", "Rprofile.site"),
#		Sys.getenv("R_PROFILE_USER"),
#		file.path(getwd(), ".Rprofile") )
#
#Filter(file.exists, candidates)