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

# Create some fuctions for later testing
f <- function(x){
  x[1]^2 + 2*x[2] + x[3]
}
# This is for a 3-dimension testing.

h <- function(x){
  x[1]^2 + 2*x[2] + x[3] + log(x[4])
}
# This is for a 4-dimension testing.

library(mvtnorm)
myNorm <- function(x){
  dmvnorm(x, mean=rep(0, dimension), sigma=diag(rep(1, dimension)))
}
# This is to compare with anthor ``correct" function.

# Write the tests
library(testthat)
library(cubature)

# Test if the result of sg.int is consistent.
test_that("The function returns a correct value.", 
          expect_equal(as.vector(sg.int(g = f, lower = c(1, 2, 3), upper = c(6, 7, 8))),
                       4557.337,
                       tolerance = 0.01))

test_that("The function returns a correct value.", 
          expect_equal(as.vector(sg.int(g = h, lower = c(1, 2, 3, 4), upper = c(5, 6, 7, 8))),
                       12146.93,
                       tolerance = 0.01))

#test_that("The answer is close to the real value.", 
#          expect_equal(as.vector(sg.int(g = myNorm, lower = c(1, 2, 3), upper = c(6, 7, 8))),
#                       adaptIntegrate(g = myNorm, lowerLimit = c(1, 2, 3), upperLimit = c(6, 7, 8))$integral,
#                       tolerance = 1))

# Test if the output is of the right class.
test_that("The output is in a matrix",
          expect_is(sg.int(g = h, lower = c(1, 2, 3, 4), upper = c(5, 6, 7, 8)), "matrix"))


# ---------- Compere the speed with and without parallel ----------

library(microbenchmark)

# compare the speed in 2 dimensions
microbenchmark(
  "Noparallel_2dim" = sg.int(g = sin, lower = c(-2, -2), upper = c(2, 2)),
  "Parallel_2dim." = sg.int(g = sin, lower = c(-2, -2), upper = c(2, 2), parallel = TRUE),
  "adaptIntegrate_2dim" = sg.int(g = sin, lower = c(-2, -2), upper = c(2, 2)),
  times = 200L
)
# No parallel is faster than parallel, and than adaptIntegrate by a small margin. Their means 
# are all around 5 to 6 milliseconds.

# compare the speed in 3 dimensions
microbenchmark(
  "Noparallel_3dim" = sg.int(g = f, lower = c(-2, -2, -2), upper = c(1, 1, 1)),
  "Parallel_3dim." = sg.int(g = f, lower = c(-2, -2, -2), upper = c(1, 1, 1), parallel = TRUE),
  "adaptIntegrate_3dim" = sg.int(g = f, lower = c(-2, -2, -2), upper = c(1, 1, 1)),
  times = 200L
)
# No parallel is faster than parallel, and than adaptIntegrate by a small margin. Their means 
# are all around 24 to 25 milliseconds.

# compare the speed in 4 dimensions
microbenchmark(
  "Noparallel_4dim" = sg.int(g = h, lower = c(-2, -2, -2, 1), upper = c(2, 3, 4, 5)),
  "Parallel_4dim." = sg.int(g = h, lower = c(-2, -2, -2, 1), upper = c(2, 3, 4, 5), parallel = TRUE),
  "adaptIntegrate_4dim" = sg.int(g = h, lower = c(-2, -2, -2, 1), upper = c(2, 3, 4, 5)),
  times = 200L
)
# No parallel is still faster than parallel by a small margin.


