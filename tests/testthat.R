shhh <- function(expr){
  capture.output(x <- suppressPackageStartupMessages(
    suppressMessages(suppressWarnings(expr))))
  invisible(x)
}


shhh(library(testthat))
shhh(library(DemoTools))
shhh(library(magrittr))
shhh(library(dplyr))
shhh(library(magrittr))
shhh(library(mockery))
shhh(library(tibble))
shhh(library(wpp2024))
test_check("DemoTools")
