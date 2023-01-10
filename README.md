# creating an objective function for the RNA folding problem
For a given ribonucleotide chain, the RNA folding problem consists in finding the native fold
among a very large number of possible conformations.
The native fold being the one of the lowest Gibbs free energy,the objective function should be ana estimator for this energy.

In this work we propose a tool to predict the Gibbs free energy of RNA structure using linear interpolation. 

#Training script
This script takes as input a list of pdb files
