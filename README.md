# CPSC312_Proj1
Gaussian Process Regression (GPR), written in Haskell, written generally, and more specifically applied to predicting the energy of Hydronium (H3O+) based on energies at specific molecular configurations. 

Hydronium data can be found here: https://scholarblogs.emory.edu/bowman/potential-energy-surfaces/

*note* only use some of the data there, as our optimization routine is very simple (just a linear system solver), so we could only work with a subset. We used 3000 test samples (18000 lines). 

*note* I was asked to link to, instead of rehost, the data. It is in a .xyz format, but the program (see gpr_h3o.hs) assumes it to be converted to .txt format and moved to the preceding folder.

In GPR we assume that the prior of any regression function is described by a Gaussian process. We leave out the details here, but for a description of the technique, see (for example), Kevin Murphy's *Machine Learning: A Probabilistic Perspective*, chapter 15. 

## 1-D GPR
gpr_simple.hs contains code where you may alter some function and add noise to it. It selects points on this function at random as training, and then does a Gaussian Process Regression fit on evenly dispersed points. Running mean_pred after produces the graph, e.g. as show in mean_pred.png. You can also visualize the priors and posteriors by running f_prior_graph and f_posterior_graph respectively. It also computes the standard deviation at each point sampled and stores them. Predictions are very reasonable and GPR works as intended. The example shown is a famous curve known as the topologist's sine curve. 

## Hydronium GPR
ParseData.hs is a module for reading data from the link above (just must convert the .xyz file into .txt). It reads each line, and then calculates the following distances: H1 to H2, H1 to H3, H1 to O, H2 to H3, H2 to O, and H3 to O. This is then fed to gpr_h3o.hs which performs GPR and makes predictions.

Note, however, that we wrote this very generally, and this code would work for any arbitrary molecule (any dimensions) written in the same format. 

Since the conformation/physics is determined only by the distance of the atoms to each other (not on their absolute coordinates), this becomes a 6-dimensional regression problem. We used 3000 (6-D) points with a squared exponential kernel, performed Gaussian Process regression, and did a simple search sample of 9000 points (this would have been better done on a grid and with Monte Carlo; we took the diagonal for simplicity). We found (zpe is the predicted zero-point energy in hartree, zpe_s is the standard deviation of zpe, and conform is the point in coordinate space (i.e. distance betweeen atoms):
```
*Main> (zpe, zpe_s, conform, x, y) <- main
*Main> zpe
1.8837485136142794e-4
*Main> zpe_s
9.712718283266998e-4
*Main> conform
(1><6)
 [ 1.6000336212676844, 1.602189102338542, 0.9745700247325499, 1.591470956087909, 0.9816657038677272, 0.9845002110474843 ]

```

This is very close to the *ab-initio* prediction of about 1.6e-4 hartree. This would be improved by doing a rigorous grid-search, Monte Carlo random sampling, and by doing steepest descent instead of LinearSolve (had there been more time). 
