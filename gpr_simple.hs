-- simple Gaussian Process regression prediction
-- this fit is only 1-dimensional, but the final fit will be much more
-- I think 6-dimensional

import Numeric.LinearAlgebra

-- initialise random number generator
-- g = getStdGen
-- TODO make seed different each function call? somehow?
seed = 7 -- for reproducibility


-- apply function f across an array
-- TODO might have to change this (or at least in ker_se) to
-- agree with hmatrix data type
-- this works with Vector type
fs f = \xs -> vector [f(x) | x <- toList xs]

-- play with these parameters (and function)
-- predict this function (with some noise)
-- need to first "flatten" data TODO
f = fs sin
n_train = 10 -- number of training points
n_test :: Num p => p
n_test = 50 -- number of testing points
s :: Fractional p => p
s = 0.00005 -- noise variance (so that we don't have perfect fit), assuming gaussian

-- get a random matrix based on your seed
-- r c are number of rows and columns
-- seed is the random seed
-- dist is the distribution (e.g. Uniform), from System.Random
-- TODO not sure when this will be needed, but thought it was helpful
randMat r c seed dist = reshape c $ randomVector seed dist (r*c)

-- get random dataset for problem
-- xset is a vector of numbers uniformly sampled from -5 to 5
-- TODO generalise this, I guess
xset = 10*(randomVector seed Uniform n_train)-5
-- add noise with mean 0, std s, Gaussian distributed
-- TODO make distribution more general?
yset = f (xset + s * randomVector seed Gaussian n_train)

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
  let aa = repmat (col (toList (a^^2))) 1 (size b)
  let bb = repmat (row (toList (b^^2))) (size a) 1
  let sqdist = aa + bb - 2 * (a `outer` b)
  return (exp(-0.5 * (1/param) * sqdist))

-- k = ker_se xset xset 0.1
-- TODO: Somehow convert ker_se to Matrix R type instead of Matrix Double
-- k = matrix 10 [1..100]
-- n_train = 10
-- let s_m_iden_n = diagl (replicate n_train s)
-- let Just ch = mbChol (mTm (k + s_m_iden_n))
-- let l = cholSolve ch (k + s_m_iden_n)

-- let x_test = linspace n_test (-5,5::Double)
-- let k_t = ker_se xset x_test 0.1
-- TODO: once conversion is done from previous todo, this will work
-- let lk = linearSolve l k_t
-- let ly = linearSolve l yset
-- let mu = (tr lk) <> ly

-- let big_k = ker_se x_test x_test

-- a = vector [1..10]
-- b = vector [3,5..21]
-- k = ker_se a b 0.1
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
