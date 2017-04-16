############################################
## PS 5625 - Applied Statistical Programming
## Problem Set 6
## Author: Luwei Ying
############################################

rm(list=ls())

# Pull up the help file of package SparseGrid
library(SparseGrid)
?SparseGrid


# ---------- Change the function to allow more dimensions and parallel processing ----------

sg.int <- function(g, ..., lower, upper, parallel=FALSE, cores=4){
  # require the packages
  require("SparseGrid")
  
  # if using parallel, register the number of cores
  require("doMC")
  if(parallel){
    registerDoMC(cores=cores)
  }
  
  # check if the lengths of lower and upper are the same
  if (length(lower) != length(upper)) stop("lower and upper must be of the sme length")
  
  # record the dimension for later use; record the lower and upper bounds for each dimension
  dimension <- length(lower)
  lower <- floor(lower)
  upper <- ceiling(upper)
  
  # chek if upper is greater than lower
  if (any(lower>upper)) stop("lower must be smaller than upper")
  
  # creat a matrix of grids which contains every possible combination of lower and upper
  gridss <- as.matrix(expand.grid(lapply(1:dimension, function(x){
    seq(lower[x], upper[x]-1, by=1)
  })))
  
  # create a list, which first item is a matrix of nodes and second item is a vector of weights
  sp.grid <- createIntegrationGrid('KPU', dimension=dimension, k=5)
  
  # define the nodes and weights
  nodes <- gridss[1,]+sp.grid$nodes
  weights <- sp.grid$weights
  for (i in 2:nrow(gridss)){
    nodes <- rbind(nodes,gridss[i,]+sp.grid$nodes)
    weights<- c(weights,sp.grid$weights)
  }
  
  # run the function g for each column and weigh the results
  gx.sp <- apply(nodes, 1, g,...)
  val.sp <- gx.sp %*% weights
  
  # return the results
  val.sp
}


# ---------- Write tests for the function using Testthat ----------

library(testthat)

# Evaluate that sg.int.hi.dim should produce a similar result to adaptIntegrate
test_that("The answer is close to the real value.", 
          expect_equal(as.vector(sg.int(g = dnorm, lower = c(-1.5, 0.5, 1), upper = c(1, 3.2, 4))),
                       adaptIntegrate(dnorm, lowerLimit = c(-1.5, 0.5, 1), upperLimit = c(1, 3.2, 4))$integral,
                       tolerance = 1e-3))
# This test currently does not work. "adaptIntegrate" returns a sigle number, while sg.int returns three.


# ---------- Compere the speed with and without parallel ----------

library(microbenchmark)

# compare the speed in 2 dimensions
microbenchmark(
  "Noparallel_2dim" = sg.int(g = dnorm, lower = c(-2, -2), upper = c(1, 1)),
  "Parallel_2dim." = sg.int(g = dnorm, lower = c(-2, -2), upper = c(1, 1), parallel = TRUE),
  times = 20
)
# No parallel is faster than parallel by a small margin.

# compare the speed in 3 dimensions
microbenchmark(
  "Noparallel_2dim" = sg.int(g = dnorm, lower = c(-2, -2, -2), upper = c(1, 1, 1)),
  "Parallel_2dim." = sg.int(g = dnorm, lower = c(-2, -2, -2), upper = c(1, 1, 1), parallel = TRUE),
  times = 20
)
# Parallel is faster than no parallel by a small margin.

# compare the speed in 4 dimensions
microbenchmark(
  "Noparallel_2dim" = sg.int(g = dnorm, lower = c(-2, -2, -2, -2), upper = c(1, 1, 1, 1)),
  "Parallel_2dim." = sg.int(g = dnorm, lower = c(-2, -2, -2, -2), upper = c(1, 1, 1, 1), parallel = TRUE),
  times = 20
)
# No parallel is faster than parallel by a small margin.