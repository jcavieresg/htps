rm(list = ls())
setwd("C:/Users/Usuario/Desktop/htps")

library(devtools)
library(Rcpp)
library(RcppArmadillo)
library(pkgKitten)

compileAttributes(verbose=TRUE) # Find and register Rcpp functions 
devtools::build() # build a package
devtools::load_all() # load all functions

