from openeye import oechem
from openeye import oeomega
import argparse
from openeye import oemolprop

def gen_conf(mol):
    oemols = []
    omegaOpts = oeomega.OEOmegaOptions(oeomega.OEOmegaSampling_FastROCS)
    omega = oeomega.OEOmega(omegaOpts)
    for enantiomer in oeomega.OEFlipper(mol.GetActive(), 6, True):
        enantiomer = oechem.OEMol(enantiomer)
        ret_code = omega.Build(enantiomer)
        if ret_code == oeomega.OEOmegaReturnCode_Success:
            halfMol = oechem.OEMol(mol, oechem.OEMCMolType_HalfFloatCartesian)
            oemols.append(halfMol)
        else:
            oechem.OEThrow.Warning("%s: %s" %
                                   (enantiomer.GetTitle(), oeomega.OEGetOmegaError(ret_code)))
    return oemols

# def filter_smiles(smis):
#     '''
#     :param smis: list of smiles
#     :return: (oe graph mols, smiles)
#     '''
#     ims = oechem.oemolistream()
#     ims.SetFormat(oechem.OEFormat_SMI)
#     smiles = "\n".join(list(smis))
#     ims.openstring(smiles)
#
#     oms = oechem.oemolostream()
#     oms.SetFormat(oechem.OEFormat_SMI)
#     oms.openstring()
#
#     filt = oemolprop.OEFilter(oemolprop.OEFilterType_BlockBuster)
#
#     goods = []
#     for i, mol in enumerate(ims.GetOEGraphMols()):
#         if filt(mol):
#             oechem.OEWriteMolecule(oms, mol)
#             goods.append(i)
#     actions = str(oms.GetString().decode("utf-8"))
#     actions = actions.split("\n")
#
#     oms.close()
#     ims.close()
#     return [smis[i] for i in goods], actions

def from_smiles(in_smiles, outfile):
    '''

    :param in_smiles: a list of smile strings
    :param outfile:
    :return:
    '''
    in_smiles = "\n".join(in_smiles)

    ims = oechem.oemolistream()
    ims.SetFormat(oechem.OEFormat_SMI)
    ims.openstring(in_smiles)

    ofs = oechem.oemolostream()
    ofs.SetFormat(oechem.OEFormat_OEB)
    filt = oemolprop.OEFilter(oemolprop.OEFilterType_BlockBuster)

    if not ofs.open(outfile):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % outfile)

    for mol in ims.GetOEMols():
        if filt(mol):
            result = gen_conf(mol)
            for res in result:
                oechem.OEWriteMolecule(ofs, res)

    ofs.close()
    return 0

def from_single_smiles(in_smiles, outfile):
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, in_smiles)

    ofs = oechem.oemolostream()
    ofs.SetFormat(oechem.OEFormat_OEB)
    if not ofs.open(outfile):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % outfile)
    filt = oemolprop.OEFilter(oemolprop.OEFilterType_BlockBuster)
    if filt(mol):
        result = gen_conf(mol)
        for res in result:
            oechem.OEWriteMolecule(ofs, res)

    ofs.close()
    return 0

def from_input_file(infile, outfile):
    ifs = oechem.oemolistream()
    if not ifs.open(infile):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % infile)

    ofs = oechem.oemolostream()
    ofs.SetFormat(oechem.OEFormat_OEB)
    if not ofs.open(outfile):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % outfile)

    for mol in ifs.GetOEMols():
        result = gen_conf(mol)
        for res in result:
            oechem.OEWriteMolecule(ofs, res)

    ofs.close()
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='input file {.csv or .smi} where csv has column SMILES and TITLE', required=True)
    parser.add_argument('-o', type=str, help='output must be .oeb', required=True)
    args = parser.parse_args()
    assert(args.o.split(".")[-1] == 'oeb')
    from_input_file(args.i, args.o)