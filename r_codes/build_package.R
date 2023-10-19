rm(list = ls())

setwd("C:/Users/Usuario/Desktop/htps")

library(devtools)
library(Rcpp)
library(RcppArmadillo)
library(pkgKitten)

compileAttributes(verbose=TRUE) # Find and register Rcpp functions 
pkgbuild::compile_dll()
devtools::load_all() # load all functions


# devtools::document() # create roxygen2 document, don't make .md by hand, let roxygen2 do it
# devtools::check() # check whether the package is OK
# devtools::build() # build a package

