#!/usr/bin/env python
#Firas Akermi

import numpy as np
import argparse
import math
import pandas as pd
from training import Read_files, CalculDistance
from plot_energy import Plot_Gibbs
def Interpolation(scores,pdb_file):
    name = pdb_file[0].split("/")[-1]
    distances = Read_files(pdb_file)
    scores = Plot_Gibbs(scores)
    x_col = 'distance'
    df = scores.sort_values(by=x_col)
    df_interp = {"AA":[],"AU":[],"AC":[],"AG":[],"UU":[],"UG":[],"UC":[],"CC":[],"CG":[],"GG":[]}
    gibbs = 0
    for col in distances:
        for i in distances[col]:
    # Find the two closest x values in the dataframe to the target x value
            x_new = i
            df['delta'] = abs(df[x_col] - x_new)
            df = df.sort_values(by='delta')
            x1, x2 = df.iloc[0][x_col], df.iloc[1][x_col]
            y1, y2 = {}, {}
            y1[col] = df.iloc[0][col]
            y2[col] = df.iloc[1][col]
            df = df.drop(columns=["delta"])
            # Interpolate the y values using a linear interpolation formula
            formula= y1[col] + (y2[col] - y1[col]) * (x_new - x1) / (x2 - x1)
            gibbs+=formula
            df_interp[col].append(formula)
    print("The Predicted Gibbs free energy for the RNA structure {0} is : {1}".format(name.split(".")[0],gibbs))
    
            

    

    
    

            
         

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Python script to Plot free energy")
    parser.add_argument("-s","--Scores_directory",help="The directory tha contains scores files",type = str, required = True)
    parser.add_argument("-i","--pdb_file",help="The path to the pdb file to predict it's energy",nargs='+',type = str, required = True)
    args = parser.parse_args()
    Interpolation(args.Scores_directory,args.pdb_file)
