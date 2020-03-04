import os
os.environ['OPENMM_CPU_THREADS'] = '1'
from mpi4py import MPI
import pandas as pd
import sys
import policy
from impress_md import interface_functions
import time
import subprocess

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def collate(path, dbase_name, target_name):
    import os, time
    path_to_watch = path
    before = {}
    while 1:
        time.sleep(10)
        for f in os.listdir(path_to_watch):
            if f not in before and os.path.exists(path +"/" + f + "/done.txt"):
                try:
                    df = pd.read_csv(path +"/" +f + "/metrics.csv")
                    before[f] = df
                except FileNotFoundError:
                    pass
        try:
            ds = pd.concat(list(before.values()))
            ds['target'] = target_name
            ds['dbase'] = dbase_name
            ds.to_csv("out.csv")
            print("TOTAL", ds.shape)
        except ValueError:
            pass



def setup_server():
    status_ = MPI.Status()

    dockPolicy = policy.MasterDockPolicy()
    mmPolicy = policy.MasterMinimizePolicy()
    # print("Master setup server.")
    ts = time.time()
    ts_start = ts

    docked_count = 0
    param_count = 0
    while True:
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status_)
        if len(data) == 1: #pipeline 1
            dtype = data.pop(0)
            if dtype == 'dock':
                r = dockPolicy(data[0])
                comm.send(r, dest=status_.Get_source(), tag=11)
            elif dtype == 'min': #pipeline 2
                r = mmPolicy(data[0])
                comm.send(r, dest=status_.Get_source(), tag=11)
            elif dtype == 'mmgbsa': #pipline 3
                pass
        else:
            pass
            # print("got some weird data", data)
        if time.time() - ts > 15:
            ts = time.time()
            # print("current counts", docked_count, param_count, time.time() - ts_start)


def worker(df, path_root, docking_only=False, receptor_file=None):
    size = comm.Get_size()

    struct = "input/"
    docker,recept = interface_functions.get_receptr()

    dock_policy, min_policy = policy.DockPolicy(), policy.MinimizePolicy()

    mols_docked = 0
    mols_minimzed = 0

    for pos in range(rank - 2, df.shape[0], size - 2):
        path = path_root + str(pos) + "/"
        try:
            smiles = df.iloc[pos,0]
            name = df.iloc[pos, 1]
            r = dock_policy(smiles)
            if r:
                # print("Rank", rank, pos, "running docking...")
                score = interface_functions.RunDocking_(smiles,struct,path, dock_obj=docker, write=True, recept=recept, name=name)
                mols_docked += 1

                if not docking_only:
                    if mols_docked % 1000 == 0:
                        comm.send(['mini', min_policy.rollout()], dest=0, tag=11)
                        r = comm.recv(source=0, tag=11)
                        min_policy.collect_rollout(r)
                    r = min_policy(smiles, score)

                    if r:
                        interface_functions.ParameterizeOE(path)
                        mscore = interface_functions.RunMinimization_(path, path)
                        mols_minimzed += 1
                        # comm.send([smiles, score, mscore], dest=0, tag=11)
                        # r = comm.recv(source=0, tag=11)
                        # print("Rank", rank, "should I run mmgbsa for 1 ns given a energy minmization result of", mscore, "?\t my model says", bool(r))
                        if r:
                            escore = interface_functions.RunMMGBSA_(path,path)
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
    df = pd.read_csv(args.smiles, sep=' ', header=None)
    num_mols = df.shape[0]
    path_root = args.path

    from shutil import copyfile

    copyfile(args.receptor_file, 'input/receptor.oeb')

    if not os.path.exists(path_root):
        os.mkdir(path_root)

    if rank == 0:
        setup_server()
    elif rank == 1:
        collate(path_root, args.dbase_name, args.target_name)
    else:
        worker(df, path_root + "/rank", docking_only=args.dock_only, receptor_file=args.receptor_file)