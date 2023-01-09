#!/usr/bin/env python
#Firas Akermi
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import os
import glob

def Plot_Gibbs(input_direct):
    '''
    Function to read scores files and to plot the scores distrubtions
    '''
    all_files = glob.glob(os.path.join(input_direct , "*.txt"))
    li = []
    for filename in all_files:
        df = pd.read_csv(filename)
        li.append(df)
        frame = pd.concat(li, axis=1)
        frame["distance"] = [i for i in range(0,21)]
    return frame
    
def PlotDistr(input_dir,output_dir):
    '''
    Function to plot the scores
    '''
    scores = Plot_Gibbs(input_dir)
    for s in scores.columns:
        if s!="distance":
            plt.figure()
            plt.plot( scores["distance"],scores[s],linestyle='-', color='b')
            plt.xlabel("Distance (Ä)")
            plt.ylabel("Score")
            plt.title("Pair {}".format(s))               
            plt.savefig("{0}{1}".format(output_dir,s))

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Python script to Plot free energy")
    parser.add_argument("-i","--Scores_directory",help="The directory of scores files",type = str, required = True)
    parser.add_argument("-o","--Output_directory",help="Output directory",type = str, required = True)
    args = parser.parse_args()
    PlotDistr(args.Scores_directory,args.Output_directory)
