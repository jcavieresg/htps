rm(list = ls())
setwd("C:/Users/Usuario/Desktop/mypackage")

#setwd("C:/Users/Usuario/Desktop/example2")

library(devtools)
library(Rcpp)
library(RcppArmadillo)
library(pkgKitten)

#RcppArmadillo.package.skeleton( "mypackage" )

compileAttributes(verbose=TRUE) # Find and register Rcpp functions 
devtools::load_all() # load all functions

# #pkgbuild::compile_dll()
# devtools::document() # create roxygen2 document, don't make .md by hand, let roxygen2 do it
# devtools::check() # check whether the package is OK
# devtools::build() # build a package

