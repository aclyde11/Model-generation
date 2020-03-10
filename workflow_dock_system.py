import sys
import pandas as pd
from impress_md import interface_functions
import os

if __name__ == '__main__':
    smiles_files = pd.read_csv(sys.argv[1], sep=' ', header=None)
    target_filoe = sys.argv[2] #twenty of these
    dbase_name = 'test'
    target_name = 'pl_pro'

    output_location = 'output_test/'
    if not os.path.exists(output_location):
        os.mkdir(output_location)

    ### OPTION 1, no need to read receptor from file system everytime

    docker, receptor = interface_functions.get_receptr(target_filoe)

    for pos in range(smiles_files.shape[0]):
        smiles = smiles_files.iloc[pos, 0]
        ligand_name = smiles_files.iloc[pos, 1]

        ## WORKFLOW 0.1 , CPU Code. doesn not run on power.
        score, res = interface_functions.RunDocking_A(smiles, None, output_location + ligand_name + '/', dbase_name, target_name,
                                        pos=pos, write=True,
                                        receptor_file=None, name=ligand_name,
                                        docking_only=True, dock_obj=docker, recept=receptor)
        interface_functions.ParameterizeOE(output_location + ligand_name + '/')


        ### THIS IS GPU CODE BELOW
        ## WORKFLOW 1
        print("Got here")
        interface_functions.RunMinimization_(output_location + ligand_name + '/', output_location + ligand_name + '/', write=True, gpu=False)

        ## Workflow 1.1 # if workflow 1 is good, we then run this.
        interface_functions.RunMMGBSA_(output_location + ligand_name + '/', output_location + ligand_name + '/', gpu=False)

        print(res)


