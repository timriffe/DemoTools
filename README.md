# DemoTools

[![Build Status](https://travis-ci.org/timriffe/DemoTools.svg?branch=master)](https://travis-ci.org/timriffe/DemoTools)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/timriffe/DemoTools?branch=master&svg=true)](https://ci.appveyor.com/project/timriffe/DemoTools)
[![codecov](https://codecov.io/gh/timriffe/DemoTools/branch/master/graph/badge.svg)](https://codecov.io/gh/timriffe/DemoTools) 
[![](https://img.shields.io/badge/devel%20version-0.12.30-yellow.svg)](https://github.com/timriffe/DemoTools)
[![issues](https://img.shields.io/github/issues-raw/timriffe/DemoTools.svg)](https://github.com/timriffe/DemoTools/issues)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)

# Tools for the evaluation, adjustment, and standardization of demographic data
Date: 2018-12-15
 
This repository contains simple functions in a package format, and is in active development. This project is commissioned by the [UN Population Division](http://www.un.org/en/development/desa/population/) and financed by the [Bill and Melinda Gates Foundation](https://www.gatesfoundation.org/) as part of the [Making Family Planning Count](http://www.un.org/en/development/desa/population/projects/making-family-planning-count/index.shtml) project. Work is also done in collaboration with Sean Fennell, and [Jose Manuel Aburto](http://findresearcher.sdu.dk/portal/en/persons/jose-manuel-aburto(34dcae96-a13a-4c4d-a941-985152180869).html), with minor contributions from [several more](https://github.com/timriffe/DemoTools/graphs/contributors)(thank you!). This work is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License ([CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/)). 

If you detect a bug or have a suggestion please notify us using the [Issues](https://github.com/timriffe/DemoTools/issues) tab on github. Even better if you fix it and make a pull request! See [CONTRIBUTING.md](https://github.com/timriffe/DemoTools/blob/master/CONTRIBUTING.md) for more tips on reporting bugs or offering patches. 

You can load the ```DemoTools``` package in R like so:
```r
# install.packages("devtools")

library(devtools)
install_github("timriffe/DemoTools")
```
(if either of the first two icons at the top of this README are red, then this might not be working at the moment. You can assume we're fixing it. If they're green, then it'll probably work.)

You can then get started like so:
```
library(DemoTools)
# interesting top-level functions include:

# for age heaping:
?Whipple
?Myers
?Bachi
?CoaleLi
?Noumbissi
?Spoorenberg
?KannistoHeap (Kannisto's old-age heaping index)
?Jdanov (Jdanov's old-age heaping index)

# other age-structure quality measures:
?ageRatioScore  # methods including "UN", "Zelnick", "Ramachandran"
?sexRatioScore
?ageSexAccuracy # methods including "UN", "Zelnick", "Ramachandran", and "Das Gupta"


# Comparison methods
?IRD (index of relative difference)
?ID (index of dissimilarity)
?survRatioError

# graduation methods
?sprague
?beers # methods including "ord" and "mod", as well as johnson option for young ages
?grabill
?splitMono
?monoCloseout
?splitOscillate # accepting e.g. spragueSimple, beersSimple as split methods

# various smoothing methods
?agesmth # including Carrier-Farrag, Arriaga, Karup-King-Newton, United Nations, Strong, Zigzag, and MAV methods
?spencer
?zelnik
?mav
?zigzag

# various lifetable evaluation and calculation functions
?ADM # and ?RDM, implementing PAS LIFIT
?LTabr # with fine control over a(x) assumptions, extrapolation, and open age groups

# interpolation
?interp (arithmetic, logarithmic, power)

# adjustment
?rescaleAgeGroups (including for cases of different age groupings)
```
Note, at the moment the Sprague family of functions can return negatives in the oldest ages, and a solution is being sought. Please stand by.

These top-level functions have implied an even larger set of simple utilities, which itself is growing fast. Presently top-level + utilities = 97 documented functions, with more in development. 

Presently all functions are in a testing phase, but the aim is to end up with a set of robust generic functions around which wrappers can be easily built for various institutional data production needs. As-is, these functions may also be useful for DIY demographers. This set of methods is a cherry-pick from legacy methods collections, including PAS, DAPPS, MPCDA, MortPack, IREDA, UN Manual X, G. Feeney Spreadsheets, formulas found in Siegel and Swanson or Shyrock and Siegel, and various (apparent) first-implementations from formulas in papers, or ad hoc DIY approximations from old pros. 

## about those icons 
Every time this repository is updated the entire code base is rebuilt on a server somewhere, and undergoes a series of checks. This happens on a Linux machine and on a Windows machine. Any warnings or errors in these builds will yield a red fail tag, and successes are green passes. Code coverage indicates what percentage of lines of code undergo formal unit testing of some kind.

