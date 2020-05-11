from openeye import oechem
from openeye import oeomega
import argparse

def main(infile, outfile):

    ifs = oechem.oemolistream()
    if not ifs.open(args.i):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % infile)

    ofs = oechem.oemolostream()
    if not ofs.open(args.b):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % outfile)

    omegaOpts = oeomega.OEOmegaOptions(oeomega.OEOmegaSampling_FastROCS)
    omega = oeomega.OEOmega(omegaOpts)

    for mol in ifs.GetOEMols():
        oechem.OEThrow.Info("Title: %s" % mol.GetTitle())

        for enantiomer in oeomega.OEFlipper(mol.GetActive(), 6, True):
            enantiomer = oechem.OEMol(enantiomer)
            ret_code = omega.Build(enantiomer)
            if ret_code == oeomega.OEOmegaReturnCode_Success:
                oechem.OEWriteMolecule(ofs, enantiomer)
            else:
                oechem.OEThrow.Warning("%s: %s" %
                                       (enantiomer.GetTitle(), oeomega.OEGetOmegaError(ret_code)))

    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='input file {.csv or .smi} where csv has column SMILES and TITLE', required=True)
    parser.add_argument('-o', type=str, help='output must be .oeb', required=True)
    args = parser.parse_args()
    assert(args.o.split(".")[-1] == 'oeb')
    main(args.i, args.o)