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

    mol = oechem.OEMol()

    ifs = oechem.oemolistream(args.i)
    ifs.SetConfTest(oechem.OEAbsCanonicalConfTest())



    iterator = iter(enumerate(ifs.GetOEMols()))
    while True:
        try:
            i, mol = next(iterator)
            print(mol.GetTitle())

            if i == 10:
                break
        except StopIteration:
            break

    print()
    ifs.close()

    ifs = oechem.oemolistream(args.i)
    ifs.SetConfTest(oechem.OEAbsCanonicalConfTest())
    for i,mol in enumerate(ifs.GetOEMols()):
        print(mol.GetTitle())
        if i== 10:
            break
