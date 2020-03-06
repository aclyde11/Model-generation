import os

from mpi4py import MPI
import pandas as pd
import sys
import policy
from impress_md import interface_functions
import time
import subprocess
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
os.environ['CUDA_VISIBLE_DEVICES'] = str(rank)

def setup_server():
    status_ = MPI.Status()
    storage = {}

    dockPolicy = policy.DockPolicy()
    mmPolicy = policy.MinimizePolicy()
    print("Master setup server.")
    ts = time.time()
    docked_count = 0
    param_count = 0
    while True:

        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status_)
        if len(data) == 1:  # pipeline 1
            res = dockPolicy(*data)
            comm.send(int(res), dest=status_.Get_source(), tag=11)
            docked_count += 1
        elif len(data) == 2:  # pipeline 2
            res = mmPolicy(*data)
            comm.send(int(res), dest=status_.Get_source(), tag=11)
            param_count += 1
        elif len(data) == 3:  # pipline 3
            res = policy.mmgbsa_ns_policy(data[0], data[1], data[2])
            comm.send(int(res), dest=status_.Get_source(), tag=11)
        else:
            print("got some weird data", data)
        if time.time() - ts > 100:
            print("current counts", docked_count, param_count)


def worker(df, path_root, dbase_name, target_name, docking_only=False, receptor_file=None):
    size = comm.Get_size()
    struct = "input/"
    docker,recept = interface_functions.get_receptr(receptor_file=receptor_file)


    for pos in range(rank, df.shape[0], size):
        pos, smiles, name = pos, df.iloc[pos, 0], df.iloc[pos, 1]
        path = path_root + str(pos) + "/"
        try:
            score, res = interface_functions.RunDocking_(smiles,struct,path, dbase_name, target_name, dock_obj=docker, write=True, recept=recept, receptor_file=receptor_file, name=name, docking_only=False)
            interface_functions.ParameterizeOE(path)
            mscore = interface_functions.RunMinimization_(path, path)
            escore = interface_functions.RunMMGBSA_(path, path)
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