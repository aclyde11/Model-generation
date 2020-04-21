# Workflow 0
from openeye import oechem
import pandas as pd
import numpy as np
from impress_md import interface_functions
import argparse
import os

world_size = int(0)
rank = int(os.environ['ALPS_APP_PE'])

def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help='input csv for smiles', required=True, type=str)
    parser.add_argument("-o", help='output file for data', required=True, type=str)
    parser.add_argument("-r", help='receptor file', required=True, type=str)
    parser.add_argument('-v', help='verbose', action='store_true')
    parser.add_argument("-n", type=int, default=1)
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
    force_flipper = True
    high_resolution = True

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
            for i, col in enumerate(columns):
                value = str(smiles_file.iloc[pos, i]).strip()
                if col.lower() != 'smiles' and 'na' not in value.lower() and len(value) > 1:
                    try:
                        oechem.OESetSDData(ligand, col, value)
                    except ValueError:
                        pass
            oechem.OEWriteMolecule(ofs, ligand)

    if ofs is not None:
        ofs.close()
