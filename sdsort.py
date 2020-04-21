from openeye import oechem
import argparse
import numpy as np

def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True)
    parser.add_argument('-o', type=str, required=True)
    parser.add_argument('-t', type=str, required=True)
    return parser.parse_args()

if __name__ == '__main__':
    args = getargs()

    ifs = oechem.oemolistream(args.i)
    ofs = oechem.oemolostream(args.o)

    mol = oechem.OEGraphMol()
    mols = []
    scores = []
    while oechem.OEReadMolecule(ifs, mol):
        mols.append(oechem.OEMol(mol))
        scores.append(float(oechem.OEGetSDData(mol, args.t)))
    ifs.close()
    sort_order = np.argsort(scores)
    print(scores)
    for i in sort_order:
        oechem.OEWriteMolecule(ofs, mols[i])
    ofs.close()