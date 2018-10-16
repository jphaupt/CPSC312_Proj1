-- simple Gaussian Process regression prediction
-- this fit is only 1-dimensional, but the final fit will be much more
-- I think 6-dimensional

import Numeric.LinearAlgebra as NLA
import Graphics.Rendering.Chart.Easy
import Graphics.Rendering.Chart.Backend.Cairo
import Data.Random.Normal
import System.Random

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
f = fs (\x -> sin(1/(5*x)))
n_train = 11 -- number of training points
n_test :: Num p => p
n_test = 50 -- number of testing points
s_noise :: Fractional p => p
s_noise = 0.00005 -- noise variance (so that we don't have perfect fit), assuming gaussian
x_start :: Fractional p => p
x_start = 0.01
x_fin :: Fractional p => p
x_fin = 1

-- get a random matrix based on your seed
-- r c are number of rows and columns
-- seed is the random seed
-- dist is the distribution (e.g. Uniform), from System.Random
-- TODO not sure when this will be needed, but thought it was helpful
randMat r c seed dist = reshape c $ randomVector seed dist (r*c)

-- get random dataset for problem
-- xset is a vector of numbers uniformly sampled from -5 to 5
-- TODO generalise this, I guess
xset = (x_fin-x_start)*(randomVector seed Uniform n_train)+x_start
-- xset = NLA.vector [-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
-- add noise with mean 0, std s, Gaussian distributed
-- TODO make distribution more general?
yset = f xset + s_noise * randomVector seed Gaussian n_train

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
Just lt = mbChol (trustSym (k + s_m_iden_n))
l = tr lt
--Just ch = mbChol (trustSym (k + s_m_iden_n))
--l = cholSolve ch (k + s_m_iden_n)

-- points we make the prediction at
--x_test = linspace n_test (-5,5::Double)
x_test = linspace n_test (x_start, x_fin)

-- computing the mean at our points
k_t = ker_se xset x_test ker_val
Just lk = linearSolve l k_t
y_matrix = col (toList yset) -- (11,1)
Just ysol = linearSolve l y_matrix
mu = (tr lk) NLA.<> ysol
{-Just lk = linearSolve l k_t
y_matrix = col (toList yset)
Just ly = linearSolve l y_matrix
mu = (tr lk) NLA.<> ly-}

-- computing the variance at our test points
k_test = ker_se x_test x_test ker_val
(k_t_x, k_t_y) = size k_test
k_test_diag = col (replicate k_t_x 1)
lk_2_sum = matrix_col_sum (lk ^^ 2)
s2 = k_test_diag - lk_2_sum
s = sqrt s2

-- prepare for plots
x_y_p = pairing (toList xset) (toList yset)
x_fxt_p = pairing (toList x_test) (toList (f (x_test)))
x_mu_p = pairing (toList x_test) (toList (flatten mu))
-- plot
mean_pred = toFile def "mean_pred.png" $ do
  layout_title .= "Mean predictions"
  plot (points "original data" x_y_p)
  plot (line "x_test against f(x_test)" [x_fxt_p])
  plot (line "x_set against average" [x_mu_p])

{-




-- draw samples from the prior
new_iden = diagl (replicate n_test 0.000001)
Just prior_ch = mbChol (trustSym (new_iden))
new_l = cholSolve prior_ch (k_test + new_iden)
rand_matr = randMat n_test 10 2342432322 Gaussian
f_prior = new_l NLA.<> rand_matr

-- TODO: clean up this code if possible (tried lambda function implementation, but failed)
prior_0_col = (toList (flatten (f_prior ¿ [0])))
prior_1_col = (toList (flatten (f_prior ¿ [1])))
prior_2_col = (toList (flatten (f_prior ¿ [2])))
prior_3_col = (toList (flatten (f_prior ¿ [3])))
prior_4_col = (toList (flatten (f_prior ¿ [4])))
prior_5_col = (toList (flatten (f_prior ¿ [5])))
prior_6_col = (toList (flatten (f_prior ¿ [6])))
prior_7_col = (toList (flatten (f_prior ¿ [7])))
prior_8_col = (toList (flatten (f_prior ¿ [8])))
prior_9_col = (toList (flatten (f_prior ¿ [9])))

prior_0 = pairing (toList x_test) prior_0_col
prior_1 = pairing (toList x_test) prior_1_col
prior_2 = pairing (toList x_test) prior_2_col
prior_3 = pairing (toList x_test) prior_3_col
prior_4 = pairing (toList x_test) prior_4_col
prior_5 = pairing (toList x_test) prior_5_col
prior_6 = pairing (toList x_test) prior_6_col
prior_7 = pairing (toList x_test) prior_7_col
prior_8 = pairing (toList x_test) prior_8_col
prior_9 = pairing (toList x_test) prior_9_col

f_prior_graph = toFile def "f_prior.png" $ do
  layout_title .= "Ten samples from the GP prior"
  plot (line "prior_0" [prior_0])
  plot (line "prior_1" [prior_1])
  plot (line "prior_2" [prior_2])
  plot (line "prior_3" [prior_3])
  plot (line "prior_4" [prior_4])
  plot (line "prior_5" [prior_5])
  plot (line "prior_6" [prior_6])
  plot (line "prior_7" [prior_7])
  plot (line "prior_8" [prior_8])
  plot (line "prior_9" [prior_9])

-- draw samples from the posterior
lk_dot = (tr' lk) NLA.<> lk
Just posterior_ch = mbChol (trustSym (new_iden - lk_dot))
new_l_post = cholSolve posterior_ch (k_test + new_iden - lk_dot)
f_posterior = mu - (new_l_post NLA.<> rand_matr)

-- TODO: clean up this code if possible (tried lambda function implementation, but failed)
posterior_0_col = (toList (flatten (f_posterior ¿ [0])))
posterior_1_col = (toList (flatten (f_posterior ¿ [1])))
posterior_2_col = (toList (flatten (f_posterior ¿ [2])))
posterior_3_col = (toList (flatten (f_posterior ¿ [3])))
posterior_4_col = (toList (flatten (f_posterior ¿ [4])))
posterior_5_col = (toList (flatten (f_posterior ¿ [5])))
posterior_6_col = (toList (flatten (f_posterior ¿ [6])))
posterior_7_col = (toList (flatten (f_posterior ¿ [7])))
posterior_8_col = (toList (flatten (f_posterior ¿ [8])))
posterior_9_col = (toList (flatten (f_posterior ¿ [9])))

posterior_0 = pairing (toList x_test) posterior_0_col
posterior_1 = pairing (toList x_test) posterior_1_col
posterior_2 = pairing (toList x_test) posterior_2_col
posterior_3 = pairing (toList x_test) posterior_3_col
posterior_4 = pairing (toList x_test) posterior_4_col
posterior_5 = pairing (toList x_test) posterior_5_col
posterior_6 = pairing (toList x_test) posterior_6_col
posterior_7 = pairing (toList x_test) posterior_7_col
posterior_8 = pairing (toList x_test) posterior_8_col
posterior_9 = pairing (toList x_test) posterior_9_col

f_posterior_graph = toFile def "f_posterior.png" $ do
  layout_title .= "Ten samples from the GP posterior"
  plot (line "posterior_0" [posterior_0])
  plot (line "posterior_1" [posterior_1])
  plot (line "posterior_2" [posterior_2])
  plot (line "posterior_3" [posterior_3])
  plot (line "posterior_4" [posterior_4])
  plot (line "posterior_5" [posterior_5])
  plot (line "posterior_6" [posterior_6])
  plot (line "posterior_7" [posterior_7])
  plot (line "posterior_8" [posterior_8])
  plot (line "posterior_9" [posterior_9])
-}
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

pairing :: [a] -> [b] -> [(a, b)]
pairing xs ys = [ (x,y) | (x,y) <- zip xs ys]
