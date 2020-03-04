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

def collate(file):
    print("loadintg data")
    ranks = {}
    for i in range(comm.Get_size() - 2):
        ranks[i + 2] = []

    df = pd.read_csv(file, sep=' ', header=None)
    print("done")
    assigner = 0
    for pos in range(df.shape[0]):
        pos, smile, name = pos, df.iloc[pos,0], df.iloc[pos,1]
        ranks[assigner + 2].append((pos, smile, name))
        assigner += 1
        assigner = assigner % (comm.Get_size() - 2)

    for k, v in ranks.items():
        print("Sending data")
        comm.send(v, dest=k, tag=11)



def setup_server():
    status_ = MPI.Status()

    dockPolicy = policy.MasterDockPolicy()
    mmPolicy = policy.MasterMinimizePolicy()
    # print("Master setup server.")
    ts = time.time()
    ts_start = ts

    docked_count = 0
    param_count = 0
    with open('out_test.csv', 'w', buffering=1) as f:
        f.write("name,smiles,Dock,Dock_U,dbase,target\n")
        while True:
            data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status_)
            for line in data:
                print("wrote")
                f.write(line)

            if time.time() - ts > 15:
                ts = time.time()
                # print("current counts", docked_count, param_count, time.time() - ts_start)


def worker(path_root, dbase_name, target_name, docking_only=False, receptor_file=None):
    size = comm.Get_size()
    data = comm.recv(tag=11)
    struct = "input/"
    docker,recept = interface_functions.get_receptr()

    dock_policy, min_policy = policy.DockPolicy(), policy.MinimizePolicy()

    mols_docked = 0
    mols_minimzed = 0

    buffer = []

    for pos in data:
        pos, smiles, name = pos
        path = path_root + str(pos) + "/"
        try:

            r = dock_policy(smiles)
            if r:
                # print("Rank", rank, pos, "running docking...")
                score, res = interface_functions.RunDocking_(smiles,struct,path, dbase_name, target_name, dock_obj=docker, write=True, recept=recept, name=name, docking_only=docking_only)
                mols_docked += 1

                if docking_only:
                    if res is not None:
                        buffer.append(res)
                    if len(buffer) > 5:
                        comm.send(buffer, dest=0, tag=11)
                        buffer = []


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

    path_root = args.path

    from shutil import copyfile, SameFileError
    try:
        copyfile(args.receptor_file, 'input/receptor.oeb')
    except SameFileError:
        pass

    if not os.path.exists(path_root):
        os.mkdir(path_root)

    if rank == 0:
        setup_server()
    elif rank == 1:
        collate(args.smiles)
    else:
        worker(path_root + "/rank", args.dbase_name, args.target_name, docking_only=args.dock_only, receptor_file=args.receptor_file)