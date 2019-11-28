[<img src="logo.png" width=100 />](https://timriffe.github.io/DemoTools/)

[<img src="logo.png" align="left" width=100 />](https://timriffe.github.io/DemoTools/)

# DemoTools

[![Build Status](https://travis-ci.org/timriffe/DemoTools.svg?branch=master)](https://travis-ci.org/timriffe/DemoTools)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/timriffe/DemoTools?branch=master&svg=true)](https://ci.appveyor.com/project/timriffe/DemoTools)
[![codecov](https://codecov.io/gh/timriffe/DemoTools/branch/master/graph/badge.svg)](https://codecov.io/gh/timriffe/DemoTools) 
[![](https://img.shields.io/badge/devel%20version-01.03.01-yellow.svg)](https://github.com/timriffe/DemoTools)
[![issues](https://img.shields.io/github/issues-raw/timriffe/DemoTools.svg)](https://github.com/timriffe/DemoTools/issues)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)

# Tools for the evaluation, adjustment, and standardization of demographic data
Date: 2019-11-28
 
`DemoTools` is an R package that contains simple functions often used in demographic analysis. It is in active development. This project is commissioned by the [UN Population Division](http://www.un.org/en/development/desa/population/) and financed by the [Bill and Melinda Gates Foundation](https://www.gatesfoundation.org/) as part of the [Making Family Planning Count](http://www.un.org/en/development/desa/population/projects/making-family-planning-count/index.shtml) project. Work is also done in collaboration with Sean Fennell, and [José Manuel Aburto](https://github.com/jmaburto), [Ilya Kashnitsky](https://ikashnitsky.github.io/), [Marius Pascariu](https://github.com/mpascariu), [Jorge Cimentada](https://github.com/cimentadaj), [Monica Alexander](https://www.monicaalexander.com/), and with minor contributions from [several more](https://github.com/timriffe/DemoTools/graphs/contributors) (thank you!). This work is licensed under the Creative Commons Attribution-ShareAlike 3.0 IGO ([CC BY-SA 3.0 IGO](https://creativecommons.org/licenses/by-sa/3.0/igo/)). 

If you detect a bug or have a suggestion please notify us using the [Issues](https://github.com/timriffe/DemoTools/issues) tab on github. Even better if you fix it and make a pull request! See [CONTRIBUTING.md](https://github.com/timriffe/DemoTools/blob/master/CONTRIBUTING.md) for more tips on reporting bugs or offering patches. 


## Getting started

If you are getting started with `DemoTools` we recommend taking a look at the tutorial articles and the examples in the package documentation. 

We'll soon add a primer here, but for now you can get started by loading the package and calling up help files, which contain worknig examples to demonstrate usage and options.

You can load the ```DemoTools``` package in R like so:
```r
# install.packages("devtools")

library(devtools)
install_github("timriffe/DemoTools")
```
(if either of the first two icons at the top of this README are red, then this might not be working at the moment. You can assume we're fixing it. If they're green, then it'll probably work.)

## Note
Sometime soon there will be an overhaul of function names. We plan to switch to snake case, with method families as the first element. This is to make naming more regular and memorable, and also to activate autocomplete in RStudio or similar.

These top-level functions have implied an even larger set of simple utilities, which itself is growing fast. Presently top-level + utilities = 121 documented functions, with more in development. 

Presently all functions are in a testing phase, but the aim is to end up with a set of robust generic functions around which wrappers can be easily built for various institutional data production needs. As-is, these functions may also be useful for DIY demographers. This set of methods is a cherry-pick from legacy methods collections, including PAS, DAPPS, MPCDA, MortPack, IREDA, UN Manual X, G. Feeney Spreadsheets, formulas found in Siegel and Swanson or Shyrock and Siegel, and various (apparent) first-implementations from formulas in papers, or ad hoc DIY approximations from old pros. 

## about those icons 
Every time this repository is updated the entire code base is rebuilt on a server somewhere, and undergoes a series of checks. This happens on a Linux machine and on a Windows machine. Any warnings or errors in these builds will yield a red fail tag, and successes are green passes. Code coverage indicates what percentage of lines of code undergo formal unit testing of some kind.

