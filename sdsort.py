from openeye import oechem
import argparse
import numpy as np
from tqdm import tqdm


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True)
    parser.add_argument('-o', type=str, required=True)
    parser.add_argument('-t', type=str, required=True)
    parser.add_argument('-u', action='store_true')
    parser.add_argument('-v', default=1, type=int, choices=[0,1,2])
    parser.add_argument('-n', default=-1, type=int)

    return parser.parse_args()


if __name__ == '__main__':
    args = getargs()

    ifs = oechem.oemolistream(args.i)
    ofs = oechem.oemolostream(args.o)

    mol = oechem.OEGraphMol()

    wrote_titles = {}

    mols = []
    scores = []
    if args.v >= 1:
        pbar = tqdm(desc='reading input {}'.format(args.i))
    while oechem.OEReadMolecule(ifs, mol):
        mols.append(oechem.OEMol(mol))
        scores.append(float(oechem.OEGetSDData(mol, args.t)))

        if args.v >= 1:
            pbar.update(1)

    if args.v >=1:
        pbar.close()
    ifs.close()

    if args.v >=1:
        print("sorting...")
    sort_order = np.argsort(scores)

    for i in tqdm(sort_order,desc=('writing out {}' if not args.u else 'writing unique molname out {}').format(args.o)):
        if args.n != -1 and args.i >= args.n:
            break

        if args.u and mols[i].GetTitle() not in wrote_titles:
            oechem.OEWriteMolecule(ofs, mols[i])
            wrote_titles[mols[i].GetTitle()] = 1
        elif not args.u:
            oechem.OEWriteMolecule(ofs, mols[i])
        elif args.u and args.v >= 2:
            print("[-u] skipping ", mols[i].GetTitle())

    ofs.close()
