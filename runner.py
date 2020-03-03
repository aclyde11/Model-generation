from mpi4py import MPI
import pandas as pd
import sys
import policy
from impress_md import interface_functions
import time
import subprocess

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def setup_server():
    status_ = MPI.Status()

    dockPolicy = policy.MasterDockPolicy()
    mmPolicy = policy.MasterMinimizePolicy()
    print("Master setup server.")
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
            print("got some weird data", data)
        if time.time() - ts > 15:
            ts = time.time()
            print("current counts", docked_count, param_count, time.time() - ts_start)


def worker(df, path_root):
    size = comm.Get_size()
    struct = "input/"
    docker,recept = interface_functions.get_receptr()

    dock_policy, min_policy = policy.DockPolicy(), policy.MinimizePolicy()

    mols_docked = 0
    mols_minimzed = 0

    for pos in range(rank - 1, df.shape[0], size - 1):
        try:
            path = path_root + str(pos)  + "/"
            smiles = df.iloc[pos,0]
            r = dock_policy(smiles)
            if r:
                print("Rank", rank, "running docking...")
                score = interface_functions.RunDocking_(smiles,struct,path, dock_obj=docker, write=True, recept=recept)
                mols_docked += 1

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
                #     if r:
                #         # print("Rank", rank, "running simulation")
                #         escore = interface_functions.RunMMGBSA_(path,path)
                #         # print("Rank", rank, "ran simulation and got", escore)
        except KeyboardInterrupt:
            exit()
        except subprocess.CalledProcessError as e:
            print("Error rank", rank, e)
        except IndexError as e:
            print("Error rank", rank, e)
        except RuntimeError as e:
            print("Error rank", rank, e)

if __name__ == '__main__':
    df = pd.read_csv(sys.argv[1], sep=' ', header=None)
    num_mols = df.shape[0]
    path_root = sys.argv[2]

    if rank == 0:
        setup_server()
    else:
        worker(df, path_root)