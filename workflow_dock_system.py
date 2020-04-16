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
    path_root = 'rank'


    output_location = 'output_test/'
    if not os.path.exists(output_location):
        os.mkdir(output_location)
    path_root = output_location + path_root
    struct = "input/"
    docker, recept = interface_functions.get_receptr(receptor_file=receptor_file)

    for pos in range(smiles_files.shape[0]):
        pos, smiles, name = pos, smiles_files.iloc[pos, 0], smiles_files.iloc[pos, 1]
        path = path_root + str(pos) + "/"
        print(pos, smiles, name, path_root, path)

        # score, res = interface_functions.RunDocking_A(smiles, struct, path, dbase_name, target_name,
        #                                               dock_obj=docker, write=False, recept=recept,
        #                                               receptor_file=receptor_file, name=name, docking_only=True)
        interface_functions.RunMinimizationGAFF(path, path)
        break
        # interface_functions.ParameterizeOE(path)

        # print("here")
        #
        # mscore = interface_functions.RunMinimization_(path, path, write=True, gpu=False)
        # if mscore < -500:
        # escore = interface_functions.RunMMGBSA_(path, path, gpu=True, niter=5000) #5ns

        # collect this result string
        # with open(path + "/metrics.csv") as f:
        #     next(f)
        #     result = next(f)
        # print(score)
        # if pos == 10:
        #     break



