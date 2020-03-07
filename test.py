import os

import pandas as pd

os.environ['OPENMM_CPU_THREADS'] = '1'
from impress_md import interface_functions
import time
import glob

df = pd.read_csv("/Users/austin/Downloads/drug_bank.smi", sep=' ', header=None)

struct = "input/"
docker, recept = interface_functions.get_receptr()

start = time.time()

files = glob.glob("/Users/austin/Downloads/good_pl/rank*")
for file in files:
    interface_functions.ParameterizeOE(file + "/")

print(time.time() - start)
