# Workflow 0
import re
from openeye import oechem
import pandas as pd
import numpy as np
from impress_md import interface_functions

def get_root_protein_name(file_name):
   return file_name.split("/")[-1].split(".")[0]

def get_smiles_col(col_names):
    return int(np.where(['smile' in s.lower() for s in col_names])[0][0])

def get_ligand_name_col(col_names):
    return int(np.where(['id' in s.lower() or 'title' in s.lower() or "name" in s.lower() for s in col_names])[0][0])

if __name__ == '__main__':
    ## Use input every run
    input_smiles_file = "/Users/austin/Box/2019-nCoV/drug-screening/HTDrugScreeningLigands/ena+db.csv"
    target_file = '/Users/austin/rec_6W02_A__DU__APR_A-201.oeb.gz'  # twenty of these
    output_poses = 'example.sdf'

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
    columns = smiles_file.columns.tolist()
    smiles_col = get_smiles_col(columns)
    name_col = get_ligand_name_col(columns)

    docker, receptor = interface_functions.get_receptor(target_file, use_hybrid=use_hybrid, high_resolution=high_resolution)

    for pos in range(49997, smiles_file.shape[0]):
        smiles = smiles_file.iloc[pos, smiles_col]
        ligand_name = smiles_file.iloc[pos, name_col]
        score, res, ligand = interface_functions.RunDocking_(smiles,
                                                dock_obj=docker,
                                                pos=pos,
                                                name=ligand_name,
                                                target_name=pdb_name,
                                                force_flipper=force_flipper)
        if res is not None:
            print(res, end='')

        if ofs is not None and ligand is not None:
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