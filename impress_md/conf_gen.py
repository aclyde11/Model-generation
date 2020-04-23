import sys
from openeye import oechem, oeomega


def FromMol(mol, use_flipper=True, num_sterocenters=12, force_flipper=False):
    """
    Generates a set of conformers as an OEMol object
    Inputs:
        mol is an OEMol
        isomers is a boolean controling whether or not the various diasteriomers of a molecule are created
        num_enantiomers is the allowable number of enantiomers. For all, set to -1
    """
    omegaOpts = oeomega.OEOmegaOptions()
    omegaOpts.SetMaxConfRange("200,800")
    omegaOpts.SetRangeIncrement(8)
    omegaOpts.SetMaxSearchTime(45)
    omega = oeomega.OEOmega(omegaOpts)

    out_conf = []

    for enantiomer in oeomega.OEFlipper(mol.GetActive(), num_sterocenters, force_flipper):
        enantiomer = oechem.OEMol(enantiomer)
        ret_code = omega.Build(enantiomer)
        if ret_code == oeomega.OEOmegaReturnCode_Success:
            out_conf.append(enantiomer)

        else:
            oechem.OEThrow.Warning("%s: %s" % (mol.GetTitle(), oeomega.OEGetOmegaError(ret_code)))

    return out_conf


def FromString(smiles, use_flipper=True, force_flipper=False, num_sterocenters=6):
    """
    Generates an set of conformers from a SMILES string
    """
    mol = oechem.OEMol()
    if not oechem.OESmilesToMol(mol, smiles):
        print("SMILES invalid for string", smiles)
        return None
    else:
        return FromMol(mol, use_flipper, num_sterocenters, force_flipper=force_flipper)
