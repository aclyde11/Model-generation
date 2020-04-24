# Workflow 0
from openeye import oechem
import pandas as pd
import numpy as np
from impress_md import interface_functions
import argparse
import os

world_size = int(0)
rank = int(os.environ['ALPS_APP_PE'])

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
    parser.add_argument("-n", type=int, default=1)
    parser.add_argument("-l", type=str, default=None, required=False)
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
    world_size = args.n

    ## Use input every run
    input_smiles_file = args.i
    target_file = args.r  # twenty of these
    basename, file_ending = ".".join(args.o.split(".")[:-1]), args.o.split(".")[-1]
    output_poses = basename + "_" + str(rank) + "." + file_ending

    ## setting don't change
    use_hybrid = True
    force_flipper = False
    high_resolution = False

    # set logging if used
    ofs = oechem.oemolostream(output_poses)


    # get root receptor Name
    pdb_name = get_root_protein_name(target_file)

    smiles_file = pd.read_csv(input_smiles_file)
    columns = smiles_file.columns.tolist()
    smiles_col = get_smiles_col(columns)
    name_col = get_ligand_name_col(columns)

    docker, receptor = interface_functions.get_receptor(target_file, use_hybrid=use_hybrid,
                                                        high_resolution=high_resolution)

    poss = np.array_split(list(range(smiles_file.shape[0])), world_size)[rank]
    for pos in poss:
        try:
            with timeout(seconds=150):
                smiles = smiles_file.iloc[pos, smiles_col]
                ligand_name = smiles_file.iloc[pos, name_col]
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


    if ofs is not None:
        ofs.close()

