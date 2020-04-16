import os

import pandas as pd
import sys
import policy
from impress_md import interface_functions
import time
import subprocess
from pymol import cmd
import re

world_size = os.environ['OMPI_COMM_WORLD_SIZE']
rank = os.environ['OMPI_COMM_WORLD_RANK']
os.environ['CUDA_VISIBLE_DEVICES'] = str(rank)


def worker(df, path_root, dbase_name, target_name, docking_only=False, receptor_file=None):
    struct = "input/"
    docker, recept = interface_functions.get_receptr(receptor_file=receptor_file)

    for pos in range(int(rank), df.shape[0], int(world_size)):
        pos, smiles, name = pos, df.iloc[pos].loc['smiles'], df.iloc[pos].loc['name']
        path = path_root + str(pos) + "/"
        try:
            score, res = interface_functions.RunDocking_A(smiles, struct, path, dbase_name, target_name,
                                                          dock_obj=docker, write=False, recept=recept,
                                                          receptor_file=receptor_file, name=name, docking_only=True)
            interface_functions.ParameterizeOE(path)
            mscore = interface_functions.RunMinimization_(path, path, write=True, gpu=True)
            print(smiles, score, mscore)
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
    parser.add_argument('--data_dump', required=True, type=str)
    parser.add_argument('-n', type=int, default=10, help='number of top moleucles to pull from data dump ')
    parser.add_argument('--dock_only', action='store_true')
    parser.add_argument('--path', type=str, required=True)
    parser.add_argument('--receptor_file', type=str, required=True)
    return parser.parse_args()


if __name__ == '__main__':
    import os

    args = get_args()

    target = re.sub("receptor.oeb", "", args.receptor_file)[:-1].split("/")[-1] + "_dock"
    print("from the receptor file provided, the target name is ", target)

    path_root = args.path

    if rank == 0:
        if not os.path.exists(path_root):
            os.mkdir(path_root)
        else:
            raise ValueError(path_root, "already exists. Please rename or delete that direcotry")

    if args.data_dump.split(".")[-1] == 'csv':
        df = pd.read_csv(args.data_dump)
    elif args.data_dump.split(".")[-1] == 'pkl':
        df = pd.read_pickle(args.data_dump)
    else:
        raise ValueError("Not sure how to parse that file yet, check data dump")

    df = df[['smiles', target]].dropna().sort_values(target)
    print("loaded up", target, "with", df.shape[0], "values in it")
    df = df.iloc[:args.n]
    print("Cutting down to ", df.shape[0], "top molecules for", target)

    worker(df, path_root + "/rank", "placeholder", target, docking_only=args.dock_only,
           receptor_file=args.receptor_file)
