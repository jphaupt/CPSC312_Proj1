-- simple Gaussian Process regression prediction 
-- this fit is only 1-dimensional, but the final fit will be much more
-- I think 6-dimensional

import Numeric.LinearAlgebra

-- apply function f across an array
-- TODO might have to change this (or at least in ker_se) to 
-- agree with hmatrix data type
-- this works with Vector type
fs f = \xs -> vector [f(x) | x <- toList xs]

-- predict this function (with some noise)
-- need to first "flatten" data TODO
-- TODO add noise, where to sample, etc.
f = fs sin

-- squared exponential kernel
-- a and b are datasets, param is the kernel parameters
-- a and b are matrices, param is a real number
-- TODO fs should be column-wise! 
-- should return a matrix (num samples in a by num samples in b)
-- TODO will have to add values to the diagonal to represent noise
ker_se a b param = 
  -- TODO
  -- have to sum across columns for n-dimensional case
  -- sqdist is calculated very differently in n-dimensions
  let sqdist = (???)
  return exp(-0.5 * (1/param) * sqdist) -- won't work because sqdist is dim (n,n)
