-- Simple Gaussian Process regression prediction
-- this fit is only 1-dimensional, but the final fit will be much more

import Numeric.LinearAlgebra as NLA
import Graphics.Rendering.Chart.Easy
import Graphics.Rendering.Chart.Backend.Cairo
import Data.Random.Normal
import System.Random

-- ***** CONSTANTS *****
-- initialise random number generator
seed = 7 -- for reproducibility

-- apply function f across an array
fs f = \xs -> vector [f(x) | x <- toList xs]

-- play with these parameters (and function)
-- predict this function (with some noise)
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

-- get random dataset for problem
-- xset is a vector of numbers uniformly sampled from -5 to 5
xset = (x_fin-x_start)*(randomVector seed Uniform n_train)+x_start
-- xset = NLA.vector [-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
-- add noise with mean 0, std s, Gaussian distributed
yset = f xset + s_noise * randomVector seed Gaussian n_train

-- sample some input points and noisy versions of the function at these points
ker_val = 0.1
k = ker_se xset xset ker_val
s_noise_iden_n = diagl (replicate n_train s_noise)
Just lt = mbChol (trustSym (k + s_noise_iden_n))
l = tr lt

-- points we make the prediction at
x_test = linspace n_test (x_start, x_fin)

-- computing the mean at our points
k_t = ker_se xset x_test ker_val
Just lk = linearSolve l k_t
y_matrix = col (toList yset) -- (11,1)
Just ysol = linearSolve l y_matrix
mu = (tr lk) NLA.<> ysol

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

-- ***** PLOTTING FUNCTIONS *****
-- plot
mean_pred = toFile def "mean_pred.png" $ do
  layout_title .= "Mean predictions"
  plot (points "original data" x_y_p)
  plot (line "x_test against f(x_test)" [x_fxt_p])
  plot (line "x_test against average" [x_mu_p])

-- draw samples from the prior
new_iden = diagl (replicate n_test 0.000001)
Just prior_ch = mbChol (trustSym (k_test + new_iden))
prior_l = tr prior_ch
prior_rand_mat = randMat n_test n_train 2342432322 Gaussian
f_prior = prior_l NLA.<> prior_rand_mat

f_prior_2plot = [toList (flatten (f_prior ¿ [i])) | i <- [0..n_train-1]]

prior_2plot = [pairing (toList x_test) (f_prior_2plot !! i) | i <- [0..n_train-1]]

-- plot *all* priors
f_prior_graph = toFile def "f_prior.png" $ do
  layout_title .= "Ten samples from the GP prior"
  plot_helper "prior" prior_2plot (n_train-1)

-- draw samples from the posterior
lk_dot = (tr' lk) NLA.<> lk
Just post_ch = mbChol (trustSym (k_test + new_iden - lk_dot))
post_l = tr post_ch
post_rand_mat = randMat n_test n_train 232432322 Gaussian
f_posterior = mu - (post_l NLA.<> post_rand_mat)

f_post_2plot = [toList (flatten (f_posterior ¿ [i])) | i <- [0..n_train-1]]

post_2plot = [pairing (toList x_test) (f_post_2plot !! i) | i <- [0..n_train-1]]

-- plot *all* posteriors
f_posterior_graph = toFile def "f_posterior.png" $ do
  layout_title .= "samples from the GP posterior"
  plot_helper "posterior" post_2plot (n_train-1)


-- ***** FUNCTIONS *****
-- squared exponential kernel
-- a and b are datasets, param is the kernel parameters
-- a and b are vectors (will be matrices for n-D case),
-- param is a real number
-- should return a matrix (num samples in a by num samples in b)
ker_se :: NLA.Vector Double -> NLA.Vector Double -> NLA.Matrix Double -> NLA.Matrix R
ker_se a b param = do
  let aa = repmat (col (toList (a^^2))) 1 (size b)
  let bb = repmat (row (toList (b^^2))) (size a) 1
  let sqdist = aa + bb - 2 * (a `outer` b)
  (exp (-0.5 * (1/param) * sqdist))

-- It will produce a matrix that sums up all the columns together
-- m is a Matrix R
matrix_col_sum :: NLA.Matrix R -> NLA.Matrix R
matrix_col_sum m = do
  let (_, c) = size m
  let f_list = map (\i -> sum_col_helper m i) [0..(c-1)]
  col f_list

-- Sums up a specific column in the matrix
sum_col_helper m i = do
  let col = m ¿ [i]
  sum (toList (flatten col))

-- pairs elements in two lists
pairing :: [a] -> [b] -> [(a, b)]
pairing xs ys = [ (x,y) | (x,y) <- zip xs ys]

-- get a random matrix based on your seed
-- r c are number of rows and columns
-- seed is the random seed
-- dist is the distribution (e.g. Uniform), from System.Random
randMat r c seed dist = reshape c $ randomVector seed dist (r*c)

-- plot first n items of lst (set of pairs like prior_2plot)
-- and call each plot str0, str1, etc.
plot_helper str lst 0 = plot (line (str ++ (show 0)) [lst !! 0])
plot_helper str lst n = do
  plot (line (str ++ (show n)) [lst !! n])
  plot_helper str lst (n-1)
