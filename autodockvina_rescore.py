
import subprocess
from openbabel.openbabel import OBConversion, OBMol

from openbabel import pybel
from tqdm import tqdm
import pickle
import multiprocessing
import argparse
from openeye import oechem

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vina', type=str, default='/Users/austin/Downloads/autodock_vina_1_1_2_mac_catalina_64bit/bin/vina')
    parser.add_argument('--receptor', type=str, required=True)
    parser.add_argument('--ligands', type=str, required=True)
    parser.add_argument('--out', type=str, required=True)
    return parser.parse_args()



def get_aff(x):
    x = str(x).split()
    next = False
    for i in x:
        if next:
            return float(i.strip())
        if "Affinity:" in i:
            next = True
    return None

def runvina(infile, outfile, tmp_file='test.pdbqt'):
    obconversion = OBConversion()
    obconversion.SetInFormat("sdf")
    obconversion.SetOutFormat("pdbqt")
    obmol = OBMol()
    notatend = obconversion.ReadFile(obmol, infile)
    obmol2 = OBMol(obmol)

    ofs = pybel.Outputfile("sdf", outfile, overwrite=True)
    pbar = tqdm()

    while notatend:
        pbar.update(1)
        if obconversion.WriteFile(obmol, 'test.pdbqt'):
            try:
                x = subprocess.check_output(["/Users/austin/Downloads/autodock_vina_1_1_2_mac_catalina_64bit/bin/vina", "--score_only", "--receptor", "/Users/austin/adrp.pdbqt" ,
                                             "--ligand", "test.pdbqt"], shell=False)
                mol2 = pybel.Molecule(obmol2)
                mol2.data.update({'AutodockVinaRescore' : str(get_aff(x))})
                ofs.write(mol2)

            except subprocess.CalledProcessError as e:
                print(e)
            except ValueError  as e:
                print(e)
        else:
            print("error writing")

        obmol = OBMol()
        notatend = obconversion.Read(obmol)
    pbar.close()


print(runvina('/Users/austin/Box/2019-nCoV/drug-screening/RELEASES/april20/ADRP-ADPR/noncovalent/adrp_adpr_recA_sorted_top100.sdf', 'out.sdf'))

