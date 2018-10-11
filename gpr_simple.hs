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
-- a and b are vectors (will be matrices for n-D case), 
-- param is a real number
-- should return a matrix (num samples in a by num samples in b)
-- TODO will have to add values to the diagonal to represent noise
-- ker_se :: Vector Double -> Vector Double -> Double -> Matrix Double
-- param = 0.1
-- a = vector [1..10]
-- b = vector [3,5..21]
ker_se a b param = do
  -- TODO will have to change this significantly for n-dimensional case
  let aa = repmat (col (toList (a^^2))) 1 10
  let bb = repmat (row (toList (b^^2))) 10 1
  let sqdist = aa + bb - 2 * (a `outer` b)
  return (exp(-0.5 * (1/param) * sqdist))

-- transpose
-- a = matrix 3 [1..9]
-- sumElements isn't what I want right now
-- need to probably implement something simliar to np.sum
-- kernel a b  = expm (sqdist)
--   where a_pow = (^^) a 2
--         b_pow = (^^) b 2
--         b_trans = (tr' b)
--         first_term = flatten (sumElements a_pow)
--         second_term = sumElements b_pow
--         third_term = (*) 2 (a Numeric.LinearAlgebra.<> b_trans)
--         sqdist = ((+) first_term second_term)
--        -- sqdist = (-) ((+) first_term second_term) third_term
--        -- exp_val = (*) ((*) -0.5 (1.0 / 2.0)) sqdist

-- http://wiki.c2.com/?DotProductInManyProgrammingLanguages
