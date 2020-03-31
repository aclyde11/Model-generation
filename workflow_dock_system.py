import sys
import pandas as pd
from impress_md import interface_functions
import os
import numpy as np

if __name__ == '__main__':
    smiles_files = pd.read_csv(sys.argv[1], sep=' ', header=None)
    receptor_file = sys.argv[2] #twenty of these
    dbase_name = 'ena_db'
    target_name = receptor_file
    path_root = 'rank/'


    output_location = 'output_test/'
    if not os.path.exists(output_location):
        os.mkdir(output_location)

    struct = "input/"
    docker, recept = interface_functions.get_receptr(receptor_file=receptor_file)

    for pos in range(smiles_files.shape[0]):
        pos, smiles, name = pos, smiles_files.iloc[pos, 0], smiles_files.iloc[pos, 1]
        path = path_root + str(pos) + "/"

        score, res = interface_functions.RunDocking_A(smiles, struct, path, dbase_name, target_name,
                                                      dock_obj=docker, write=True, recept=recept,
                                                      receptor_file=receptor_file, name=name, docking_only=False)

        interface_functions.ParameterizeOE(path)

        mscore = interface_functions.RunMinimization_(path, path, write=True, gpu=True)
        if mscore < -500:
            escore = interface_functions.RunMMGBSA_(path, path, gpu=True, niters=5000) #5ns

        #collect this result string
        with open(path + "/metrics.csv") as f:
            next(f)
            result = next(f)



