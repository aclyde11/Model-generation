# Workflow 0
from openeye import oechem
import pandas as pd
import numpy as np
from impress_md import interface_functions
import argparse
import os
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
world_size = comm.Get_size()

import signal

class timeout:
    def __init__(self, seconds=1, error_message='Timeout'):
        self.seconds = seconds
        self.error_message = error_message
    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)
    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)
    def __exit__(self, type, value, traceback):
        signal.alarm(0)

def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help='input csv for smiles', required=True, type=str)
    parser.add_argument("-o", help='output file for data', required=True, type=str)
    parser.add_argument("-r", help='receptor file', required=True, type=str)
    parser.add_argument('-v', help='verbose', action='store_true')
    return parser.parse_args()


def get_root_protein_name(file_name):
    return file_name.split("/")[-1].split(".")[0]


def get_smiles_col(col_names):
    return int(np.where(['smile' in s.lower() for s in col_names])[0][0])


def get_ligand_name_col(col_names):
    return int(np.where(['id' in s.lower() or 'title' in s.lower() or "name" in s.lower() for s in col_names])[0][0])


def dockStructure(data):
    pos, smiles, ligand_name, pdb_name, force_flipper = data
    score, res, ligand = interface_functions.RunDocking_(smiles,
                                                         dock_obj=docker,
                                                         pos=pos,
                                                         name=ligand_name,
                                                         target_name=pdb_name,
                                                         force_flipper=force_flipper)
    # print(oechem.OEWriteMolToBytes(".smi", ligand))

    return score, res, ligand


if __name__ == '__main__':
    args = getargs()
    input_smiles_file = args.i
    target_file = args.r  # twenty of these
    basename, file_ending = ".".join(args.o.split(".")[:-1]), args.o.split(".")[-1]
    output_poses = basename + "_" + str(rank) + "." + file_ending
    pdb_name = get_root_protein_name(target_file)
    ## setting don't change
    use_hybrid = True
    force_flipper = False
    high_resolution = True
    ofs = oechem.oemolostream(output_poses)

    if rank == 0:
        smiles_file = pd.read_csv(input_smiles_file)
        columns = smiles_file.columns.tolist()
        smiles_col = get_smiles_col(columns)
        name_col = get_ligand_name_col(columns)

        for pos in range(0, smiles_file.shape[0], 5):
            poss = list(range(pos, pos+5))
            smiles = smiles_file.iloc[pos : min(pos+5, smiles_file.shape[0]), smiles_col]
            ligand_name = smiles_file.iloc[pos : min(pos+5,smiles_file.shape[0]), name_col]
            status = MPI.Status()
            comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            rank_from = status.Get_source()
            data = zip(poss, smiles, ligand_name)
            comm.send(data, dest=rank_from, tag=23)




    else:
        docker, receptor = interface_functions.get_receptor(target_file, use_hybrid=use_hybrid,
                                                            high_resolution=high_resolution)

        comm.send([], dest=0, tag=11)
        poss = comm.recv(source=0, tag=MPI.ANY_TAG)

        while True:
            for (pos, smiles, ligand_name) in poss:
                try:
                    with timeout(seconds=150):
                        score, res, ligand = interface_functions.RunDocking_(smiles,
                                                                     dock_obj=docker,
                                                                     pos=pos,
                                                                     name=ligand_name,
                                                                     target_name=pdb_name,
                                                                     force_flipper=force_flipper)

                        if args.v:
                            print("RANK {}:".format(rank), res, end='')
                        if ofs and ligand is not None:
                            oechem.OEWriteMolecule(ofs, ligand)
                except TimeoutError:
                    print("TIMEOUT", smiles, ligand_name)
                    continue
            comm.send([], dest=0, tag=11)
            poss = comm.recv(source=0, tag=MPI.ANY_TAG)


        if ofs is not None:
            ofs.close()

