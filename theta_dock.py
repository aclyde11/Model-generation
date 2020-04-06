# Workflow 0

import sys
import pandas as pd
from impress_md import interface_functions
from openeye import oechem

if __name__ == '__main__':
    use_conformer_file = True
    target_filoe = sys.argv[2] #twenty of these
    dbase_name = 'test'
    target_name = 'pl_pro'

    if use_conformer_file:
        docker, receptor = interface_functions.get_receptr(target_filoe, has_ligand=False)

        ifs = oechem.oemolistream("/Users/austin/Downloads/out.oeb")
        ifs.SetConfTest(oechem.OEOmegaConfTest())

        mols = []
        for mol in ifs.GetOEMols():
            mol = oechem.OEMol(mol)
            mols.append(mol)

        for pos, mol in enumerate(mols):
            smiles = oechem.OECreateSmiString(mol)
            ligand_name = mol.GetTitle()

            score = interface_functions.RunDocking_preconf(mol, dock_obj=docker)

            res = "{},{},{},{},{},{},{}\n".format(str(pos), ligand_name, smiles, score, 0, dbase_name, target_name)
            print(res)

    else:
        smiles_files = pd.read_csv(sys.argv[1], sep=' ', header=None)
        target_filoe = sys.argv[2] #twenty of these
        dbase_name = 'test'
        target_name = 'pl_pro'

        docker, receptor = interface_functions.get_receptr(target_filoe)

        for pos in range(smiles_files.shape[0]):
            smiles = smiles_files.iloc[pos, 0]
            ligand_name = smiles_files.iloc[pos, 1]
            score = interface_functions.RunDocking_(smiles, dock_obj=docker)

            res = "{},{},{},{},{},{},{}\n".format(str(pos), ligand_name, smiles, score, 0, dbase_name, target_name)
            print(res)