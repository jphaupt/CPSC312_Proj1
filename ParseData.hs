module ParseData where

import Data.List.Split
import Numeric.LinearAlgebra
import Data.String.Utils

calculateLines s lines =
  do
    let y = convertToFloat (head (splitWs (lines !! (s+1))))
    let h1 = tail (splitWs (lines !! (s+2)))
    let h2 = tail (splitWs (lines !! (s+3)))
    let h3 = tail (splitWs (lines !! (s+4)))
    let o = tail (splitWs (lines !! (s+5)))
    return [y, calDist h1 h2, calDist h1 h3, calDist h1 o, calDist h2 h3, calDist h2 o, calDist h3 o]


-- Given a file, it will parse the data and return a matrix
-- The resulting matrix will be (total_lines/6)*7
-- The first column will be the y values
-- The following six will be the Euclidean distance of the following
--  d(H1, H2), d(H1, H3), d(H1, O), d(H2, H3), d(H2, O), d(H3,O)
readcsv filename =
  do
    file <- readFile filename
    let lines = splitOn "\n" file
    let indexes = [0, 0+6 .. ((length lines)-2)]
    let r = [ atom |  i <- indexes, atom <- calculateLines i lines]
    return (fromLists r)

-- try:
-- readcsv "test.csv"

convertToFloat s = read s :: Float
calDist [x1, y1, z1] [x2, y2, z2] = sqrt (x'*x' + y'*y' + z'*z')
    where
      x' = convertToFloat(x1) - convertToFloat(x2)
      y' = convertToFloat(y1) - convertToFloat(y2)
      z' = convertToFloat(z1) - convertToFloat(z2)
