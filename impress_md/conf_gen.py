import sys
from openeye import oechem, oeomega

def FromMol(mol, isomer=True, num_enantiomers=1):
    """
    Generates a set of conformers as an OEMol object
    Inputs:
        mol is an OEMol
        isomers is a boolean controling whether or not the various diasteriomers of a molecule are created
        num_enantiomers is the allowable number of enantiomers. For all, set to -1
    """
    omegaOpts = oeomega.OEOmegaOptions()
    omegaOpts.SetMaxConfs(25)
    omega = oeomega.OEOmega(omegaOpts)
    out_conf = []
    
    if not isomer:
        ret_code = omega.Build(mol)
        if ret_code == oeomega.OEOmegaReturnCode_Success:
            out_conf.append(mol)
        else:
            oechem.OEThrow.Warning("%s: %s" % (mol.GetTitle(), oeomega.OEGetOmegaError(ret_code)))
    
    elif isomer:
        for enantiomer in oeomega.OEFlipper(mol.GetActive(),6,True):
            enantiomer = oechem.OEMol(enantiomer)
            ret_code = omega.Build(enantiomer)
            if ret_code == oeomega.OEOmegaReturnCode_Success:
                out_conf.append(enantiomer)
                num_enantiomers -= 1
                if num_enantiomers == 0:
                    break
            else:
                oechem.OEThrow.Warning("%s: %s" % (mol.GetTitle(), oeomega.OEGetOmegaError(ret_code)))

    return out_conf

def FromString(smiles, isomer=True, num_enantiomers=1):
    """
    Generates an set of conformers from a SMILES string
    """
    mol = oechem.OEMol()
    if not oechem.OESmilesToMol(mol,smiles):
        print("SMILES invalid for string", smiles)
        return None
    else:
        return FromMol(mol,isomer,num_enantiomers)

# TODO: Placeholder
def SelectEnantiomer(mol_list):
    return mol_list[0]
