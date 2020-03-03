import pandas as pd
import sys
import policy
from impress_md import interface_functions
import time
import subprocess

df = pd.read_csv("input/john_smiles_kinasei.smi", sep=' ')

struct = "input/"
docker, recept = interface_functions.get_receptr()

start = time.time()
for pos in range(0, 5):
    path = "test" + str(pos)  + "/"
    smiles = df.iloc[pos,0]
    score = interface_functions.RunDocking_(smiles,struct,path, dock_obj=docker, recept=recept)
    interface_functions.ParameterizeOE(path)
    mscore = interface_functions.RunMinimization_(path, path)
print(time.time() - start)