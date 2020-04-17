from impress_md import interface_functions
import subprocess
import re
import pandas as pd
import os
import shutil
from glob import glob



def worker(df, path_root, dbase_name, target_name, docking_only=False, receptor_file=None, oe=False):
    struct = "input/"

    for pos in range(int(rank), df.shape[0], int(world_size)):
        pos, smiles, name = pos, df.iloc[pos].loc['smiles'], "index_" + str(df.index.tolist()[pos])
        path = path_root + str(pos) + "/"
        try:
            docker, recept = interface_functions.get_receptr(receptor_file=receptor_file)

            score, res = interface_functions.RunDocking_A(smiles, struct, path, dbase_name, target_name,
                                                          dock_obj=docker, write=False, recept=recept,
                                                          receptor_file=receptor_file, name=name, docking_only=True, oe=oe)
            if not docking_only:
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
    parser.add_argument('--param', action='store_true')
    parser.add_argument('--path', type=str, required=True)
    parser.add_argument('--receptor_file', type=str, required=True)
    parser.add_argument('--agg_file', type=str, default=None)
    parser.add_argument('--overwrite', action='store_true', help='allow overwrite of dir')
    args = parser.parse_args()
    args.agg = not(args.agg_file is None)
    print(args)
    return args



if __name__ == '__main__':
    rank = int(os.environ['OMPI_COMM_WORLD_RANK'])

    args = get_args()

    if args.agg:
        from openeye import oechem, oedocking
        from mpi4py import MPI
        if rank == 0:
            print("Using agg requires mpi4py.")
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        world_size = comm.Get_size()
    else:
        comm = None
        world_size = os.environ['OMPI_COMM_WORLD_SIZE']
    os.environ['CUDA_VISIBLE_DEVICES'] = str(rank)

    target = "ADRP-ADPR_pocket1_round1_dock"

    if rank == 0:
        print("from the receptor file provided, the target name is ", target)

    path_root = args.path
    if rank == 0:
        if not os.path.exists(path_root):
            print("made path root", path_root)
            os.mkdir(path_root)
        elif args.overwrite:
            print("deleting and overwriting", path_root)
            shutil.rmtree(path_root, ignore_errors=True)
            os.mkdir(path_root)
        else:
            raise ValueError(path_root, "already exists. Please rename or delete that direcotry")

    if args.data_dump.split(".")[-1] == 'csv':
        df = pd.read_csv(args.data_dump)
    elif args.data_dump.split(".")[-1] == 'pkl':
        df = pd.read_pickle(args.data_dump)
    else:
        raise ValueError("Not sure how to parse that file yet, check data dump")

    # dfcols = df.columns.tolist()
    # dfcols = list(map(lambda x : re.sub('_round1', '', x), dfcols))
    # print(dfcols)
    # df.columns = dfcols
    df = df[['smiles', target]].dropna().sort_values(target)
    if rank == 0:
        print("loaded up", target, "with", df.shape[0], "values in it")
    df = df.iloc[:args.n]

    if rank == 0:
        print("Cutting down to ", df.shape[0], "top molecules for", target)

    worker(df, path_root + "/rank", "placeholder", target, docking_only=(not args.param),
           receptor_file=args.receptor_file, oe=args.agg)

    if args.agg:
        comm.Barrier()
        if rank == 0:
            print("now aggregating...")

            files = glob(path_root + "/*/lig.oeb")

            protein = oechem.OEGraphMol()
            oedocking.OEReadReceptorFile(protein, args.receptor_file)
            omolout = oechem.oemolostream(args.agg_file + "com_out.pdb")
            oechem.OEWriteMolecule(omolout, protein)
            omolout.close()

            omolout = oechem.oemolostream(args.agg_file)
            omolin = oechem.oemolistream()
            for file in files:
                lig = oechem.OEGraphMol()
                omolin.open(file)
                oechem.OEReadMolecule(omolin, lig)
                oechem.OEWriteMolecule(omolout, lig)
            omolout.close()
            omolin.close()
        comm.Barrier()
