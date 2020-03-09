import sys
import pandas as pd
from impress_md import interface_functions


if __name__ == '__main__':
    smiles_files = pd.read_csv(sys.argv[1], sep=' ', header=None)
    target_filoe = sys.argv[2]
    dbase_name = 'test'
    target_name = 'pl_pro'

    WHICH_OPTION=0
    ### OPTION 1, no need to read receptor from file system everytime

    if WHICH_OPTION:
        docker, receptor = interface_functions.get_receptr(target_filoe)

        for pos in range(smiles_files.shape[0]):
            smiles = smiles_files.iloc[pos, 0]
            ligand_name = smiles_files.iloc[pos, 1]
            score, res = interface_functions.RunDocking_(smiles, None, None, dbase_name, target_name,
                                            pos=pos, write=True,
                                            receptor_file=None, name=ligand_name,
                                            docking_only=True, dock_obj=docker, recept=receptor)
            print(res)


    ### OPTION 2, need to read receptor eevery time
    else:

        for pos in range(smiles_files.shape[0]):
            smiles = smiles_files.iloc[pos, 0]
            ligand_name = smiles_files.iloc[pos, 1]
            score, res = interface_functions.RunDocking_(smiles, target_filoe, None, dbase_name, target_name,
                                                         pos=pos, write=True,
                                                         receptor_file=target_filoe, name=ligand_name,
                                                         docking_only=True, dock_obj=None, recept=None)
            print(res)