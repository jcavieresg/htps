rm(list = ls())
setwd("C:/Users/Usuario/Desktop/htps")

library(devtools)
library(Rcpp)
library(RcppArmadillo)
library(pkgKitten)

compileAttributes(verbose=TRUE) # Find and register Rcpp functions 
pkgbuild::compile_dll() # compile the functions
devtools::load_all() # load all functions

