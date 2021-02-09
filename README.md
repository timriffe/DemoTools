[<img src="logo.png" align="left" width=100 />](https://timriffe.github.io/DemoTools/)

# DemoTools

[![R build status](https://github.com/timriffe/DemoTools/workflows/R-CMD-check/badge.svg)](https://github.com/timriffe/DemoTools/actions)
[![codecov](https://codecov.io/gh/timriffe/DemoTools/branch/master/graph/badge.svg)](https://codecov.io/gh/timriffe/DemoTools) 
[![](https://img.shields.io/badge/devel%20version-01.13.08-yellow.svg)](https://github.com/timriffe/DemoTools)
[![issues](https://img.shields.io/github/issues-raw/timriffe/DemoTools.svg)](https://github.com/timriffe/DemoTools/issues)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

# Tools for aggregate demographic analysis
Date: 2021-02-09
 
`DemoTools` is an R package that contains simple functions often used in demographic analysis. It is in active development. 

This project is commissioned by the [UN Population Division](http://www.un.org/en/development/desa/population/) and financed by the [Bill and Melinda Gates Foundation](https://www.gatesfoundation.org/) as part of the [Making Family Planning Count](http://www.un.org/en/development/desa/population/projects/making-family-planning-count/index.shtml) project. Work is also done in collaboration with Sean Fennell, [José Manuel Aburto](https://github.com/jmaburto), [Ilya Kashnitsky](https://ikashnitsky.github.io/), [Marius Pascariu](https://github.com/mpascariu), [Jorge Cimentada](https://github.com/cimentadaj), [Monica Alexander](https://www.monicaalexander.com/), and with minor contributions from [several more](https://github.com/timriffe/DemoTools/graphs/contributors) (thank you!). This work is licensed under the Creative Commons Attribution-ShareAlike 3.0 IGO ([CC BY-SA 3.0 IGO](https://creativecommons.org/licenses/by-sa/3.0/igo/)). 

The idea behind `DemoTools` is to provide  a common set of functions that can be easily used by analysts and scientists working on demographic analysis and modelling. 

If you detect a bug or have a suggestion please notify us using the [Issues](https://github.com/timriffe/DemoTools/issues) tab on github. Even better if you fix it and make a pull request! See [CONTRIBUTING.md](https://github.com/timriffe/DemoTools/blob/master/CONTRIBUTING.md) for more tips on reporting bugs or offering patches. 


## <i class="fa fa-cog" aria-hidden="true"></i> Getting Started and installation

If you are getting started with `DemoTools` we recommend taking a look at the tutorial articles and the examples in the package documentation. 


You can load the ```DemoTools``` package in R like so:
```r
# install.packages("devtools")

library(devtools)
# requires the development version of rstan, sorry!
install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
install_github("timriffe/DemoTools")
```

## <i class="fa fa-pencil" aria-hidden="true"></i> Citation

To cite `DemoTools` in publications, please use:

Riffe T, Aburto JM, Alexander M,Fennell S, Kashnitsky I, Pascariu M and Gerland P. (2019) DemoTools: An R package of tools for aggregate demographic analysis
   URL: https://github.com/timriffe/DemoTools/. 
  
A BibTeX entry for LaTeX users is:

```
@Misc{demotools,
  Title		= {Demo{T}ools: {A}n {R} package of tools for aggregate demographic analysis},
  Author    = {Riffe, T and Aburto, JM and Alexander, M and Fennell, S and Kashnitsky, I and Pascariu, M and Gerland, P},
  Year      = {2019},
  note		= {URL:~\url{https://github.com/timriffe/DemoTools/}}
}

```

## <i class="fa fa-arrow-alt-circle-up" aria-hidden="true"></i> About top icons
If either of the first two icons at the top of this README are red, then the installation might not be working. You can assume we're fixing it. If they're green, then it should work.

Every time this repository is updated the entire code base is rebuilt on a server somewhere, and undergoes a series of checks. This happens on a Linux machine and on a Windows machine. Any warnings or errors in these builds will yield a red fail tag, and successes are green passes. Code coverage indicates what percentage of lines of code undergo formal unit testing of some kind.

