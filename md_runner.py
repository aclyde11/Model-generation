import os

import pandas as pd
import sys
import policy
from impress_md import interface_functions
import time
import subprocess
from pymol import cmd

world_size = os.environ['OMPI_COMM_WORLD_SIZE']
rank = os.environ['OMPI_COMM_WORLD_RANK']
os.environ['CUDA_VISIBLE_DEVICES'] = str(rank)



def worker(df, path_root, dbase_name, target_name, docking_only=False, receptor_file=None):
    struct = "input/"
    docker,recept = interface_functions.get_receptr(receptor_file=receptor_file)


    for pos in range(rank, df.shape[0], world_size):
        pos, smiles, name = pos, df.iloc[pos, 0], df.iloc[pos, 1]
        path = path_root + str(pos) + "/"
        try:
            score, res = interface_functions.RunDocking_A(smiles, struct, path, dbase_name, target_name,
                                                          dock_obj=docker, write=False, recept=recept,
                                                          receptor_file=receptor_file, name=name, docking_only=True)
            interface_functions.ParameterizeOE(path)
            mscore = interface_functions.RunMinimization_(path, path, write=True, gpu=True)
            print(smiles, score, mscore, escore)
            if mscore < -500:
                escore = interface_functions.RunMMGBSA_(path, path, gpu=True, niter=5000)  # 5n
                print(smiles, score, mscore, escore)
        except KeyboardInterrupt:
            exit()
        except subprocess.CalledProcessError as e:
            print("Error rank", rank, e)
        except IndexError as e:
            print("Error rank", rank, e)
        except RuntimeError as e:
            print("Error rank", rank, e)
        with open(path + "done.txt", 'w') as f:
            f.write("t")

def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--dock_only', action='store_true')
    parser.add_argument('--smiles', type=str, required=True)
    parser.add_argument('--path', type=str, required=True)
    parser.add_argument('--target_name', type=str, required=True)
    parser.add_argument('--dbase_name', type=str, required=True)
    parser.add_argument('--receptor_file', type=str, required=True)
    return parser.parse_args()

if __name__ == '__main__':
    import os
    args = get_args()

    path_root = args.path
    df = pd.read_csv(args.smiles, sep=' ', header=None)

    if rank == 0:
        if not os.path.exists(path_root):
            os.mkdir(path_root)
    worker(df, path_root + "/rank", args.dbase_name, args.target_name, docking_only=args.dock_only, receptor_file=args.receptor_file)
