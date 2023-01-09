#!/usr/bin/env python
#Firas Akermi

import numpy as np
import argparse
import math
import pandas as pd
from training import Read_files, CalculDistance
from plot_energy import Plot_Gibbs
def Interpolation(scores,pdb_file,pred_output):
    scores = Plot_Gibbs(scores)
    output = "{0}{1}.pred".format(pred_output,pdb_file.split("/")[-1].split(".")[0])
    print(output)
    print(scores)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Python script to Plot free energy")
    parser.add_argument("-s","--Scores_directory",help="The directory tha contains scores files",type = str, required = True)
    parser.add_argument("-i","--pdb_file",help="The path to the pdb file to predict it's energy",type = str, required = True)
    parser.add_argument("-o","--Output_directory",help="Output directory for predicted energies",type = str, required = True)
    args = parser.parse_args()
    Interpolation(args.Scores_directory,args.pdb_file,args.Output_directory)
