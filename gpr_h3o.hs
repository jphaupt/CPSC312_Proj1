-- Gaussian Process regression prediction for 6-D
-- Only using a subset of H3O data

import Numeric.LinearAlgebra as NLA
import Graphics.Rendering.Chart.Easy
import Graphics.Rendering.Chart.Backend.Cairo
import Data.Random.Normal
import System.Random
import ParseData

-- ***** CONSTANTS *****
-- initialise random number generator
seed = 7 -- for reproducibility

-- apply function f across an array
fs f = \xs -> vector [f(x) | x <- toList xs]

-- play with these parameters (and function)
-- predict this function (with some noise)
f = fs (\x -> sin(1/(5*x)))
n_test :: Num p => p
n_test = 50 -- number of testing points
s_noise :: Fractional p => p
s_noise = 0.00005 -- noise variance (so that we don't have perfect fit), assuming gaussian
ker_val = 0.1

-- ***** MAIN FUNCTION *****
main = do
  -- read data from Krems group (NOTE: small subset!)
  dat <- readcsv "../h3o_snippet.txt"
  let d = d1-1 where (_, d1) = size dat
  let xset = dat ¿ [1..d]
  let yset = dat ¿ [0]
  let (n_train, _) = size xset

  -- sample some input points and noisy versions of the function at these points
  let k = ker_se xset xset ker_val

  let s_m_iden_n = diagl (replicate n_train s_noise)
  let Just lt = mbChol (trustSym (k + s_m_iden_n))
  let l = tr lt

  let x_maxs = [maximum (toList (flatten (xset ¿ [i]))) | i <- [0..d-1]]
  let x_mins = [minimum (toList (flatten (xset ¿ [i]))) | i <- [0..d-1]]

  -- points we make the prediction at
  let x_test = fromColumns [linspace n_test (x_maxs !! i, x_mins !! i) | i <- [0..d-1]]

  -- computing the mean at our points
  let k_t = ker_se xset x_test ker_val
  let Just lk = linearSolve l k_t
  let Just ysol = linearSolve l yset
  let mu = (tr lk) NLA.<> ysol

  -- computing the variance at our test points
  let k_test = ker_se x_test x_test ker_val
  let lk_2_sum = matrix_col_sum (lk ^^ 2)
  let s2 = asColumn (takeDiag k_test) - lk_2_sum
  let s = sqrt s2
  let s_l = toList (flatten s)

  let mu_l = toList (flatten mu)
  let zpe = minimum mu_l
  -- . is function composition
  let arg_zpe = head $ filter ((== zpe) . (mu_l !!)) [0..]

  -- returns predicted zero-point energy, its standard deviation, the conformation,
  -- and (to compare) the conformations of all input with their energies
  return (zpe, s_l !! arg_zpe, xset ? [arg_zpe], xset, yset)

-- ***** HELPER FUNCTIONS *****
-- squared exponential kernel
-- a and b are datasets, param is the kernel parameters
-- a and b are matrices,
-- param is a real number
-- returns a matrix (num samples in a by num samples in b)
ker_se a b param = do
  let a_sum = flatten (matrix_row_sum (a^^2))
  let b_sum = flatten (matrix_row_sum (b^^2))
  let aa = repmat (col (toList a_sum)) 1 (size b_sum)
  let bb = repmat (row (toList b_sum)) (size a_sum) 1
  let sqdist = aa + bb - 2 * (a NLA.<> (tr b))
  (exp (-0.5 * (1/param) * sqdist))

-- produces a matrix that sums up all the columns together
-- m is a Matrix R
matrix_col_sum :: NLA.Matrix R -> NLA.Matrix R
matrix_col_sum m = do
  let (_, c) = size m
  let f_list = map (\i -> sum_col_helper m i) [0..(c-1)]
  col f_list

-- produces a matrix that sums up all the row together
matrix_row_sum :: NLA.Matrix R -> NLA.Matrix R
matrix_row_sum m = do
  let (r, _) = size m
  let f_list = map (\i -> sum_row_helper m i) [0..(r-1)]
  col f_list

-- Sums up a specific row in the matrix
sum_row_helper m i = do
  let col = m ? [i]
  sum (toList (flatten col))

-- Sums up a specific column in the matrix
sum_col_helper m i = do
  let col = m ¿ [i]
  sum (toList (flatten col))

pairing :: [a] -> [b] -> [(a, b)]
pairing xs ys = [ (x,y) | (x,y) <- zip xs ys]

-- get a random matrix based on your seed
-- r c are number of rows and columns
-- seed is the random seed
-- dist is the distribution (e.g. Uniform), from System.Random
randMat r c seed dist = reshape c $ randomVector seed dist (r*c)
