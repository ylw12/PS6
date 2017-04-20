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
  if (length(lower) != length(upper)) stop("lower and upper must be of the same length")
  
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
  x[1] + x[2]^2 + x[3]^3
}
# This is for a 3-dimension testing.

h <- function(x){
  x[1] + x[2]^2 + x[3]^3 + x[4]^4 + x[5]^5
}
# This is for a 5-dimension testing.

# Write the tests
library(testthat)
library(cubature)

# Test if the result of sg.int is correct for 3 dimensions.
test_that("The function returns a correct value.", 
          expect_equal(as.vector(sg.int(g = f, lower = c(0, 0, 0), upper = c(3, 3, 3))),
                       1215/4,
                       tolerance = 0.01))
# The function f and h we defined all can be intergrated by hand. I put the latex code 
# for the intergrating procedure of our first example. The other examples are similar:
#\begin{align*}
#&\int_0^3 \int_0^3 \int_0^3 (x_1 +x^2 +x^3) dx_1 dx_2 dx_3\\
#=&\int_0^3 \int_0^3 [(x_1x_3 +x^2 +x^3) | _0 ^3] dx_1 dx_2\\
#=&\int_0^3 [(3x_1x_2 +x_2^3 + \frac{81}{4}x_2) | _0 ^3] dx_1\\
#=&\frac{81}{2}+\frac{1053}{4}\\
#=&\frac{1215}{4}
#\end{align*}

# Test if the result of sg.int is correct for 3 dimensions when the vectors of 
# upper bound and lower bound contain different numbers.
test_that("The function returns a correct value.", 
          expect_equal(as.vector(sg.int(g = f, lower = c(0, 1, 2), upper = c(6, 7, 8))),
                       28915.29,
                       tolerance = 0.01))

# Test if the result of sg.int is correct for 5 dimensions in parallel.
test_that("The function returns a correct value.", 
          expect_equal(as.vector(sg.int(g = h, lower = c(-1, -1, -1, -1, -1), upper = c(1, 1, 1, 1, 1), parallel=TRUE)),
                       256/15,
                       tolerance = 0.01))

# Test if the input is reasonable.
expect_error(sg.int(g = f, lower = c(0, 0) , upper = c(3, 4, 5)),
             'lower and upper must be of the same length')

# Test if the output is of the right class.
test_that("The output is in a matrix",
          expect_is(sg.int(g = h, lower = c(1, 2, 3, 4), upper = c(5, 6, 7, 8)), "matrix"))


# ---------- Compere the speed of different methods ----------

library(microbenchmark)

# compare the speed in 2 dimensions
microbenchmark(
  "Noparallel_2dim" = sg.int(g = sin, lower = c(-2, -2), upper = c(2, 2)),
  "Parallel_2dim." = sg.int(g = sin, lower = c(-2, -2), upper = c(2, 2), parallel = TRUE),
  "adaptIntegrate_2dim" = adaptIntegrate(sin, lowerLimit = c(-2, -2), upperLimit = c(2, 2)),
  times = 200L
)
# No parallel is faster than parallel by a small margin. adaptIntegrate is much faster.

# compare the speed in 3 dimensions
microbenchmark(
  "Noparallel_3dim" = sg.int(g = f, lower = c(-2, -2, -2), upper = c(1, 1, 1)),
  "Parallel_3dim." = sg.int(g = f, lower = c(-2, -2, -2), upper = c(1, 1, 1), parallel = TRUE),
  "adaptIntegrate_3dim" = adaptIntegrate(f, lowerLimit = c(-2, -2, -2), upperLimit = c(1, 1, 1)),
  times = 200L
)
# No parallel is faster than parallel by a small margin. adaptIntegrate is still much faster.

# compare the speed in 5 dimensions
microbenchmark(
  "Noparallel_5dim" = sg.int(g = h, lower = c(-1, -1, -1, -1, -1), upper = c(1, 1, 1, 1, 1)),
  "Parallel_5dim." = sg.int(g = h, lower = c(-1, -1, -1, -1, -1), upper = c(1, 1, 1, 1, 1), parallel = TRUE),
  "adaptIntegrate_5dim" = adaptIntegrate(h, lowerLimit = c(-1, -1, -1, -1, -1), upperLimit = c(1, 1, 1, 1, 1)),
  times = 200L
)
# This time, using parallel with sg.int is faster than not using parallel.
# adaptIntegrate is still the fastest.


# ---------- Compere the accurancy of different methods ----------
sg.int(g = f, lower = c(0, 0, 0), upper = c(3, 3, 3)) - 1215/4
adaptIntegrate(f, lowerLimit = c(0, 0, 0), upperLimit = c(3, 3, 3))$integral - 1215/4
# adaptIntegrate is more accurant.


# ---------- Find Maximization of Functions ----------

# We want to find the local maxima for this function:

Tri <- function(x, y){
  sin(((x^2)/2) - ((y)/4))*cos(2*x - exp(y))
}

# There are (more than) two methods to do this. First try optimize(). Since our function 
# has two dimension, we are going to intergrate along one dimension first and then
# move to the other dimension.
Int_y <- function(y){
  optimize(Tri, y=y, 
           lower = c(-1,-1), upper = c(3,3), 
           maximum = T)
}

# Go along the x dimension from -1 to 3.
Allvalues <- lapply(seq(-1,3,0.0001), Int_y)

# find the value of x corresponding to the local maximum.
Max_x <- max(unlist(lapply(Allvalues, function(x) x$maximum)))

# record the index of that maximum.
Max <- which(unlist(lapply(Slices, function(x) x$maximum)) == Max_x)

# Find the corresponding Max_y.
Max_y <- seq(-1,3,0.0001)[Max]

# Find the corresponding Maximum value.
Slices[Max][[1]]$objective

# The function finds an x of 2.03, y of 1.43 and the maximum value is 0.98.

# Seond, try optim(). Modify the function and combine x and y to a single vector.
Tri <- function(x){
  sin(((x[1]^2)/2) - ((x[2])/4))*cos(2*x[1] - exp(x[2]))
}

# Optimize function over both parameters
optim(par = c(2, 1.6), fn = Tri,
      lower = c(-1,-1), upper = c(3,3), method = 'L-BFGS-B',
      control=list(fnscale=-1))

# The function finds an x of 1.95, y of 1.36 and the maximum value is 1.

# Change the initial value
optim(par = c(0, 1), fn = Tri,
      lower = c(-1,-1), upper = c(3,3), method = 'L-BFGS-B',
      control=list(fnscale=-1))

# The estimation becomes far less accurant probabily because our initial value of x is
# too far from the real value.