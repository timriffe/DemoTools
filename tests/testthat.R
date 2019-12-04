shhh <- function(expr){
  capture.output(x <- suppressPackageStartupMessages(
    suppressMessages(suppressWarnings(expr))))
  invisible(x)
}


shhh(library(testthat))
shhh(library(DemoTools))
shhh(library(magrittr))
shhh(library(dplyr))

test_check("DemoTools")
