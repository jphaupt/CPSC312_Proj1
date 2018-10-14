-- simple Gaussian Process regression prediction
-- this fit is only 1-dimensional, but the final fit will be much more
-- I think 6-dimensional

import Numeric.LinearAlgebra as NLA
import Graphics.Rendering.Chart.Easy
import Graphics.Rendering.Chart.Backend.Cairo

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
s_noise :: Fractional p => p
s_noise = 0.00005 -- noise variance (so that we don't have perfect fit), assuming gaussian

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
yset = f (xset + s_noise * randomVector seed Gaussian n_train)

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
-- ker_se :: Monad m => Vector R -> Vector R -> Matrix R -> m (Matrix R)
ker_se :: NLA.Vector Double -> NLA.Vector Double -> NLA.Matrix Double -> NLA.Matrix R
ker_se a b param = do
  -- TODO will have to change this significantly for n-dimensional case
  let aa = repmat (col (toList (a^^2))) 1 (size b)
  let bb = repmat (row (toList (b^^2))) (size a) 1
  let sqdist = aa + bb - 2 * (a `outer` b)
  (exp (-0.5 * (1/param) * sqdist))

-- sample some input points and noisy versions of the function at these points
ker_val = 0.1
k = ker_se xset xset ker_val
s_m_iden_n = diagl (replicate n_train s_noise)
Just ch = mbChol (mTm (k + s_m_iden_n))
l = cholSolve ch (k + s_m_iden_n)

-- points we make the prediction at
x_test = linspace n_test (-5,5::Double)

-- computing the mean at our points
k_t = ker_se xset x_test 0.1
Just lk = linearSolve l k_t
y_matrix = col (toList yset)
Just ly = linearSolve l y_matrix
mu = (tr lk) NLA.<> ly

-- computing the variance at our test points
k_test = ker_se x_test x_test ker_val
(k_t_x, k_t_y) = size k_test
k_test_diag = col (replicate k_t_x 1)
lk_2_sum = matrix_col_sum (lk ^^ 2)
s2 = k_test_diag - lk_2_sum
s = sqrt s2

-- Helper functions
-- It will produce a matrix that sums up all the columns together
-- m is a Matrix R
matrix_col_sum :: NLA.Matrix R -> NLA.Matrix R
matrix_col_sum m = do
  let (_, c) = size m
  let f_list = map (\i -> sum_helper m i) [0..(c-1)]
  col f_list

-- Sums up a specific column in the matrix
sum_helper m i = do
  let col = m ¿ [i]
  sum (toList (flatten col))

-- TODO: organize

pairing :: [a] -> [b] -> [(a, b)]
pairing xs ys = [ (x,y) | (x,y) <- zip xs ys]

x_y_p = pairing (toList xset) (toList yset)
x_fxt_p = pairing (toList xset) (toList (f (x_test)))
x_mu_p = pairing (toList xset) (toList (flatten mu))

mean_pred = toFile def "mean_pred.png" $ do
  layout_title .= "Mean predictions"
  plot (points "original data" x_y_p)
  plot (points "x_set against x_test" x_fxt_p)
  plot (points "x_set against average" x_mu_p)


new_iden = diagl (replicate n_test 0.000001)
Just prior_ch = mbChol (mTm (k_test + new_iden))
new_l = cholSolve prior_ch (k_test + new_iden)
-- f_prior =

-- prior = toFile def "prior.png" $ do
--   layout_title .= "Prior"
--   plot (points "original data" x_y_p)
--   plot (points "x_set against x_test" x_fxt_p)
--   plot (points "x_set against average" x_mu_p)
