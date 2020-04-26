
import subprocess
from openbabel.openbabel import OBConversion, OBMol
import numpy as np
from openbabel import pybel
from tqdm import tqdm
import pickle
import multiprocessing
import argparse
from openeye import oechem

# def get_args():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--vina', type=str, default='/Users/austin/Downloads/autodock_vina_1_1_2_mac_catalina_64bit/bin/vina')
#     parser.add_argument('--receptor', type=str, required=True)
#     parser.add_argument('--ligands', type=str, required=True)
#     parser.add_argument('--out', type=str, required=True)
#     return parser.parse_args()



def get_aff(x):
    x = str(x).split()
    next = False
    for i in x:
        if next:
            return float(i.strip())
        if "Affinity:" in i:
            next = True
    return None

def runvina(infile, outfile, receptor, tmp_file='test.pdbqt'):
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
        if obconversion.WriteFile(obmol, tmp_file):
            try:
                x = subprocess.check_output(["/Users/austin/Downloads/autodock_vina_1_1_2_mac_catalina_64bit/bin/vina", "--score_only", "--receptor", receptor ,
                                             "--ligand", tmp_file], shell=False)
                mol2 = pybel.Molecule(obmol2)
                mol2.data.update({'AutodockVinaRescoreOnly' : str(get_aff(x))})
                ofs.write(mol2)

            except subprocess.CalledProcessError as e:
                print(e)
                ofs.write(obmol)

            except ValueError  as e:
                print(e)
                ofs.write(obmol)

        else:
            print("error writing")

        obmol = OBMol()
        notatend = obconversion.Read(obmol)
        obmol2 = OBMol(obmol)

    pbar.close()
    print("FAILED")

def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', action='store_true')
    parser.add_argument('-i', default='/Users/austin/Box/2019-nCoV/drug-screening/raw_data/V5_docking_data_apri_24/NSP15/3/NSP15_3_6w01_cat_sorted.sdf')
    parser.add_argument('-o', default='/Users/austin/Box/2019-nCoV/drug-screening/raw_data/V5_docking_data_apri_24/NSP15/3/NSP15_3_6w01_cat_sorted_top100_rescore.sdf')
    parser.add_argument('-r', default="/Users/austin/Box/2019-nCoV/drug-screening/receptorsV5/NSP15_6w01_apo.pdbqt")
    parser.add_argument('-w', default=14)
    parser.add_argument('-t', default='tmp')
    return parser.parse_args()

if __name__ == "__main__":
    args = get_args()
    if args.a:
        runvina(args.i, args.o, args.r, tmp_file=args.t)
    else:
        import subprocess

        ifn = args.i
        ofn = args.o
        receptor = args.r
        tmpdir = args.t
        workers = args.w

        mols = []

        ifs = oechem.oemolistream(ifn)
        mol_ = oechem.OEMol()
        while oechem.OEReadMolecule(ifs, mol_):
            mols.append(oechem.OEMol(mol_))
            if len(mols) > 10000:
                break
        ifs.close()

        splits = np.array_split(list(range(len(mols))), workers)

        procs = []
        for i, arr in enumerate(splits):
            ofs = oechem.oemolostream(tmpdir + "/input_" + str(i) + ".sdf")
            for j in arr:
                oechem.OEWriteMolecule(ofs, mols[j])
            ofs.close()

            procs.append(subprocess.Popen(["python", "autodockvina_rescore.py", "-r", args.r, '-a', '-i',tmpdir + "/input_" + str(i) + ".sdf", '-o', tmpdir + "/output_" + str(i) + ".sdf", '-t', tmpdir + "/tmp_" + str(i) + ".pdbqt" ]))
        for proc in procs:
            proc.communicate()
