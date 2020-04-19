# Workflow 0
import re
from openeye import oechem
import pandas as pd
import numpy as np
from impress_md import interface_functions

def get_root_protein_name(file_name):
   return file_name.split("/")[-1]

def get_smiles_col(col_names):
    return int(np.where(['smile' in s.lower() for s in col_names])[0][0])

def get_ligand_name_col(col_names):
    return int(np.where(['id' in s.lower() or 'title' in s.lower() or "name" in s.lower() for s in col_names])[0][0])

if __name__ == '__main__':
    ## Use input every run
    input_smiles_file = "/Users/austin/high.csv"
    target_file = '/Users/austin/rec_6W02_A__DU__APR_A-201.oeb.gz'  # twenty of these
    output_poses = 'example.csv'

    ## setting don't change
    use_hybrid = True
    force_flipper = True
    high_resolution = True

    #set logging if used
    if output_poses is not None:
        ofs = oechem.oemolostream(output_poses)
    else:
        ofs = None

    # get root receptor Name
    pdb_name = get_root_protein_name(target_file)

    smiles_file = pd.read_csv(input_smiles_file)
    smiles_col = get_smiles_col(smiles_file.columns.tolist())
    name_col = get_ligand_name_col(smiles_file.columns.tolist())

    docker, receptor = interface_functions.get_receptor(target_file, use_hybrid=use_hybrid, high_resolution=high_resolution)

    for pos in range(smiles_file.shape[0]):
        smiles = smiles_file.iloc[pos, smiles_col]
        ligand_name = smiles_file.iloc[pos, name_col]
        score, res, ligand = interface_functions.RunDocking_(smiles,
                                                dock_obj=docker,
                                                pos=pos,
                                                name=ligand_name,
                                                target_name=pdb_name,
                                                force_flipper=force_flipper)
        if ofs is not None:
            oechem.OEWriteMolecule(ofs, ligand)


    if ofs is not None:
        ofs.close()