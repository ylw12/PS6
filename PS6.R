############################################
## PS 5625 - Applied Statistical Programming
## Problem Set 6
## Author: Luwei Ying
############################################

rm(list=ls())

# Pull up the help file of package SparseGrid
library(SparseGrid)
?SparseGrid

# ---------- Change the function to allow more dimensions ----------

sg.int <- function(g, ..., lower, upper){
  # require the package
  require("SparseGrid")
  
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


# ---------- Test the function ----------

sg.int(g=dnorm, lower = c(-1.5, 0.5, 1), upper = c(1, 3.2, 4))
