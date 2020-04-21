import argparse
from openeye import oechem
from tqdm import tqdm


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, type=str)
    parser.add_argument("-o", required=True, type=str)
    return parser.parse_args()


if __name__ == '__main__':
    args = getargs()

    ifs = oechem.oemolistream(args.i)
    ofs = oechem.oemolostream(args.o)
    lig = oechem.OEGraphMol()

    pbar = tqdm()
    while oechem.OEReadMolecule(ifs, lig):
        oechem.OEWriteMolecule(ofs, lig)
        pbar.update(1)

    pbar.close()
    ifs.close()
    ofs.close()
