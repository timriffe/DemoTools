# DemoTools

[![Build Status](https://travis-ci.org/timriffe/DemoTools.svg?branch=master)](https://travis-ci.org/timriffe/DemoTools)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/timriffe/DemoTools?branch=master&svg=true)](https://ci.appveyor.com/project/timriffe/DemoTools)
[![codecov](https://codecov.io/gh/timriffe/DemoTools/branch/master/graph/badge.svg)](https://codecov.io/gh/timriffe/DemoTools) 

Tools for the evaluation, adjustment, and standardization of demographic data

This repository contains simple functions in a package format, and is in active development. This project is commissioned by the [UN Population Division](http://www.un.org/en/development/desa/population/) and financed by the [Bill and Melinda Gates Foundation](https://www.gatesfoundation.org/) as part of the [Making Family Planning Count](http://www.un.org/en/development/desa/population/projects/making-family-planning-count/index.shtml) project. Work is also done in collaboration with Sean Fennell of the US Census Bureau. This work is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License ([CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/)). If you detect a bug or have a suggestion please notify us using the Issues tab on github. Even better if you fix it and make a pull request!

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
?Whipple
?Myers
?Bachi
?CoaleLi
?Noumbissi
?Spoorenberg
?ageRatioScore  # methods including "UN", "Zelnick", "Ramachandran"
?sexRatioScore
?ageSexAccuracy # methods including "UN", "Zelnick", "Ramachandran", and "Das Gupta"
?AHI
?WI
?IRD
?ID
?survRatioError
?spragueSimple
?grabill
?splitMono
?monoCloseout
?spragueOscillate
?beersModSimple
?birthCohorts
?adjustAge
?ADM
?RDM
# last updated 4-Aug-2017
```
These top-level functions have implied an even larger set of simple utilities, which itself is growing fast. We won't bother demonstrating those. 

Presently all functions are in first drafts, but the aim is to end up with a set of robust generic functions around which wrappers can be easily built for various institutional data production needs. As-is, these functions may also be useful for DIY demographers. This set of methods is a cherry-pick from legacy methods collections, including PAS, DAPPS, MPCDA, MortPack, IREDA, UN Manual X, G. Feeney Spreadsheets, formulas found in Siegel and Swanson or Shyrock and Siegel, and various (apparent) first-implementations from formulas in papers, or ad hoc DIY approximations from old pros. 

## about those icons 
Every time this repository is updated the entire code base is rebuilt on a server somewhere, and undergoes a series of checks. This happens on a linux machine and on a windows machine. Any warnings or errors in these builds will yield a red fail tag, and successes are green passes. Code coverage indicates what percentage of lines of code undergo testing of some kind. We've not yet set up units tests, but when we start to do so, the aim will be to push this percentage as high as possible.