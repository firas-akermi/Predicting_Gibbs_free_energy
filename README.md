# Creation of an objective function for the RNA folding problem
The RNA folding problem is the process of determining the specific three-dimensional structure of a ribonucleotide chain, called the "native fold," from the vast number of possible conformations. The native fold is the conformation with the lowest Gibbs free energy, and the goal is to find an estimator for this energy to determine the native fold.

# Training script
The training script is designed to perform the following tasks:
- Compute interatomic distances between C3' atoms from the given dataset of PDB files.
- Generate 10 distance distributions for the 10 base pairs (AA, AU, AC, AG, UU, UC, UG, CC, CG, GG), considering only distances within a single chain.
- Consider only residues separated by at least 3 positions on the sequence (i.e., residues i and i+4, i and i+5, etc.).
- Compute observed frequencies, which are 10 × 20 distances intervals (from 0 to 20 Å)
- Compute the reference frequency for the "XX" pair
- Compute the log-ratio of the two frequencies; the observed and reference frequencies.
# Plot_energy script
The plot_energy script script  will plot the interaction profiles (with R, ggplot, etc.): the score as a function of the distance.
# Prediction script
It  calculates the distances for a given RNA structure. Using the same distance thresholds of 20 Å and residues i, i+4. For each distance, the script computex a score value using linear interpolation. The estimated Gibbs free energy of the RNA conformation  is calculated by summing up all the scores from the distances.

# Dependencies
To run the scripts numpy and pandas must be installed:
```
pip install pandas
pip install numpy
```
# Pipeline execution
1. First, clone the github repository:
```
git clone https://github.com/firas-akermi/Predicting_Gibbs_free_energy.git
```
2. Make the scripts executable:
```
chmod a+x training.py plot_energy.py prediction.py
```
3. Run the training script:
```
./training -i <pdb files> -o <output directory>
```
The training script creates 10 files; Each files contains the scores of a given pair.

4. Visualizing the distances as function of scores:
```
./plot_energy.py -i <Scores directory> -o <Output directory>
```
5. Predicting Gibbs free energy:
```
./prediction.py -s <Scores directory> -i <pdb file of the RNA structure to predict>
```