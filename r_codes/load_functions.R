rm(list = ls())
setwd("C:~")

library(devtools)
library(Rcpp)
library(RcppArmadillo)
library(pkgKitten)

#RcppArmadillo.package.skeleton( "mypackage" )

compileAttributes(verbose=TRUE) # Find and register Rcpp functions 
devtools::load_all() # load all functions