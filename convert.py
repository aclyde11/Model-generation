import argparse
from openeye import oechem
from tqdm import tqdm


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, type=str)
    parser.add_argument("-o", required=True, type=str)
    parser.add_argument('-u', action='store_true')
    parser.add_argument('-v', default=1, type=int, choices=[0,1,2])
    return parser.parse_args()


if __name__ == '__main__':
    args = getargs()

    ifs = oechem.oemolistream(args.i)
    ofs = oechem.oemolostream(args.o)
    lig = oechem.OEGraphMol()
    wrote_titles = {}

    if args.v >= 1:
        pbar = tqdm(desc='reading {} and writing {}'.format(args.i, args.o))
    while oechem.OEReadMolecule(ifs, lig):
        if args.u and lig.GetTitle() not in wrote_titles:
            oechem.OEWriteMolecule(ofs, lig)
            wrote_titles[lig.GetTitle()] = 1
        elif not args.u:
            oechem.OEWriteMolecule(ofs, lig)
        elif args.u and args.v >= 2:
            print("[-u] skipping ", lig.GetTitle())
        if args.v >= 1:
            pbar.update(1)

    if args.v >= 1:
        pbar.close()
    ifs.close()
    ofs.close()
