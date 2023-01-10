#!/usr/bin/env python
#Firas Akermi
import numpy as np
import argparse
import math
import pandas as pd
def Read_files(file_name):
    '''
    function to read pdb file and to calculate distance between carbon alpha atoms
    '''
    Atoms = []
    for files in file_name:
        for line in open(files,'r'):
            if line.startswith("ATOM"):
                Atoms.append(line)
    dist = {"AA":[],"AU":[],"AC":[],"AG":[],"UU":[],"UG":[],"UC":[],"CC":[],"CG":[],"GG":[]}
    for elt in range(len(Atoms)):
        for r in range(elt,len(Atoms)):
            if abs(float(Atoms[elt][22:26].strip())-float(Atoms[r][22:26].strip())) > 3  \
            and Atoms[elt][13:16].strip() == "C3'" and Atoms[r][13:16].strip() == "C3'" :
                a = Atoms[elt][18:20].strip()+Atoms[r][18:20].strip()
                this_dist = CalculDistance(x1=float(Atoms[elt][31:38].strip()),x2 = float(Atoms[r][31:38].strip()),
                                            y1= float(Atoms[elt][39:46].strip()),y2= float(Atoms[r][39:46].strip()),
                                            z1=float(Atoms[elt][47:54].strip()),z2=float(Atoms[r][47:54].strip()))
                if this_dist<21:
                    if a =="AA":
                        dist["AA"].append(this_dist)
                    elif a == "AU" or a == "UA":
                        dist["AU"].append(this_dist)
                    elif a == "AC" or a == "CA":
                        dist["AC"].append(this_dist)
                    elif a == "AG" or a == "GA":
                        dist["AG"].append(this_dist)
                    elif a == "UU":
                        dist["UU"].append(this_dist)
                    elif a == "UG" or a == "GU":
                        dist["UG"].append(this_dist)
                    elif a == "UC" or a == "CU":
                        dist["UC"].append(this_dist)
                    elif a == "CC":
                        dist["CC"].append(this_dist)
                    elif a == "CG" or a == "GC":
                        dist["CG"].append(this_dist)
                    elif a == "GG":
                        dist["GG"].append(this_dist)
                    else:
                        raise Exception("Are you sure this is an RNA structure?")
                    
    return dist 
   
def CalculDistance(x1,x2,y1,y2,z1,z2):
    '''
    function to calculate euclidian distance
    '''
    euclidian =math.sqrt((x1 - x2)**2 + (y1-y2)**2+(z1- z2)**2)
                            
    return euclidian    

def CalculFrequencies(input_file):
    
    '''
    Function to calculate observed frequencies, refrence frequencies and scores
    '''
    distance = Read_files(input_file)
    Nr = {
        "AA":{"0":0,"1":0,"2":0,"3":0,"4":0,"5":0,"6":0,"7":0,"8":0,"9":0,"10":0,"11":0,"12":0,"13":0,"14":0,"15":0,"16":0,"17":0,"18":0,"19":0,"20":0},
        "AU":{"0":0,"1":0,"2":0,"3":0,"4":0,"5":0,"6":0,"7":0,"8":0,"9":0,"10":0,"11":0,"12":0,"13":0,"14":0,"15":0,"16":0,"17":0,"18":0,"19":0,"20":0},
        "AC":{"0":0,"1":0,"2":0,"3":0,"4":0,"5":0,"6":0,"7":0,"8":0,"9":0,"10":0,"11":0,"12":0,"13":0,"14":0,"15":0,"16":0,"17":0,"18":0,"19":0,"20":0},
        "AG":{"0":0,"1":0,"2":0,"3":0,"4":0,"5":0,"6":0,"7":0,"8":0,"9":0,"10":0,"11":0,"12":0,"13":0,"14":0,"15":0,"16":0,"17":0,"18":0,"19":0,"20":0},
        "UU":{"0":0,"1":0,"2":0,"3":0,"4":0,"5":0,"6":0,"7":0,"8":0,"9":0,"10":0,"11":0,"12":0,"13":0,"14":0,"15":0,"16":0,"17":0,"18":0,"19":0,"20":0},
        "UG":{"0":0,"1":0,"2":0,"3":0,"4":0,"5":0,"6":0,"7":0,"8":0,"9":0,"10":0,"11":0,"12":0,"13":0,"14":0,"15":0,"16":0,"17":0,"18":0,"19":0,"20":0},
        "UC":{"0":0,"1":0,"2":0,"3":0,"4":0,"5":0,"6":0,"7":0,"8":0,"9":0,"10":0,"11":0,"12":0,"13":0,"14":0,"15":0,"16":0,"17":0,"18":0,"19":0,"20":0},
        "CC":{"0":0,"1":0,"2":0,"3":0,"4":0,"5":0,"6":0,"7":0,"8":0,"9":0,"10":0,"11":0,"12":0,"13":0,"14":0,"15":0,"16":0,"17":0,"18":0,"19":0,"20":0},
        "CG":{"0":0,"1":0,"2":0,"3":0,"4":0,"5":0,"6":0,"7":0,"8":0,"9":0,"10":0,"11":0,"12":0,"13":0,"14":0,"15":0,"16":0,"17":0,"18":0,"19":0,"20":0},
        "GG":{"0":0,"1":0,"2":0,"3":0,"4":0,"5":0,"6":0,"7":0,"8":0,"9":0,"10":0,"11":0,"12":0,"13":0,"14":0,"15":0,"16":0,"17":0,"18":0,"19":0,"20":0},
    }
    N = pd.DataFrame.from_dict(Nr)
    for interval in range(0,21):
        for key1,key2 in zip(distance,N.columns):
                if key1 == key2:
                    for v in distance[key1]:
                        if v >= interval and v<= interval +1:
                            N.iloc[interval][key2]+=1
    Nxx = [i for i in N.sum(axis=0)]
    Nij = [i for i in N.sum(axis=1)]
    obs_freq = pd.DataFrame.from_dict(Nr)
    Ref_freq = pd.DataFrame.from_dict(Nr)
    score = pd.DataFrame.from_dict(Nr)
    Ref_freq=N.divide(Nxx,axis=1)
    obs_freq=N.divide(Nij, axis=0)
    for freq  in obs_freq:
            score[freq]= -1 * (np.log10(obs_freq[freq]/Ref_freq[freq]))
    for s in score.columns:        
        score.loc[score[s] >= 10, s] = 10
    score.fillna(10,inplace=True)
    for s in score.columns:
            score[s].to_csv(r'./data/scores/{}.txt'.format(s),index=False,header = True)
    
       
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="training script")
    parser.add_argument("-i","--PDB_file_path",help="input PDB file",nargs='+',type = str, required = True)
    args = parser.parse_args()
    CalculFrequencies(args.PDB_file_path)


