# CPSC312_Proj1
Gaussian Process Regression, written in Haskell, written generally, and more specifically applied to predicting the energy of Hydronium (H3O+) based on energies at specific molecular configurations. 

## Goals (TODO): 
 - fit GPR on general, generated R^1 functions (e.g. f = sin(X)), with noise added. call this e.g. gpr_simple.hs
 - parse hydronium data from Prof Krems's group (include on github iff allowed). parse_dat.hs
 - implement a modified GPR program (probably with a more complicated optimization routine) to run on hydronium data, and predict different conformation energies, e.g. ground state. gpr_main.hs
