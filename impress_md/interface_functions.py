import os
from contextlib import contextmanager

import numpy as np
from openeye import oechem, oedocking
from . import conf_gen
from . import dock_conf


@contextmanager
def working_directory(directory):
    owd = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)



def get_receptor(receptor_file=None, use_hybrid=True, high_resolution=True):
    """

    :param receptor_file: file to .oeb or oeb.gz with prepared protein receptor
    :param use_hybrid: whether or not to override using hybrid docking method
    :param high_resolution: set resolution for search
    :return: docking objector from OE, the receptor oemol
    """
    from . import dock_conf
    from openeye import oedocking

    receptor = dock_conf.PrepareReceptorFromBinary(receptor_file)

    if oedocking.OEReceptorHasBoundLigand(receptor) and use_hybrid:
        dock_method = oedocking.OEDockMethod_Hybrid
    else:
        dock_method = oedocking.OEDockMethod_Chemgauss4

    if high_resolution:
        reso = oedocking.OESearchResolution_High
    else:
        reso = oedocking.OESearchResolution_Default

    dock = oedocking.OEDock(dock_method, reso)
    dock.Initialize(receptor)
    return dock, receptor


def CanSmi(mol, isomeric, kekule):
    """
    Returns the cannonical smile from the OEMol provided
    :param mol: OEMolBase object
    :param isomeric: force isometric
    :param kekule: use kekule cleaning
    :return: string of OESmiles
    """
    oechem.OEFindRingAtomsAndBonds(mol)
    oechem.OEAssignAromaticFlags(mol, oechem.OEAroModel_OpenEye)
    smiflag = oechem.OESMILESFlag_Canonical
    if isomeric:
        smiflag |= oechem.OESMILESFlag_ISOMERIC

    if kekule:
        for bond in mol.GetBonds(oechem.OEIsAromaticBond()):
            bond.SetIntType(5)
        oechem.OECanonicalOrderAtoms(mol)
        oechem.OECanonicalOrderBonds(mol)
        oechem.OEClearAromaticFlags(mol)
        oechem.OEKekulize(mol)

    smi = oechem.OECreateSmiString(mol, smiflag)
    return smi

def RunDocking_conf(ligand: oechem.OEMol, dock_obj: oedocking.OEDock, pos: int = None, name: str = None, target_name: str = None, force_flipper: bool = True) -> object:
    """

    :param smiles: a str representing the molecule (SMILES)
    :param dock_obj: The OEDock prepared receptor
    :param pos: int for print string
    :param name: name of ligand for printing string
    :param target_name: name of target from docking
    :param force_flipper: whether or not to flip entamiers when generating conformers
    :return: score, string result, ligand
    """
    if not dock_obj.IsInitialized():
        assert(False)

    scores = []
    ligands = []



    lig = oechem.OEMol()
    dock_conf.DockConf_(dock_obj, ligand, lig, MAX_POSES=1)
    scores.append(dock_obj.ScoreLigand(lig))
    ligands.append(lig)

    # get best score from the Enantiomers
    if len(scores) > 0:
        bs = np.argmin(scores)
        score_min = scores[bs]
        ligand_min = ligands[bs]

        dockMethod = dock_obj.GetName()
        oedocking.OESetSDScore(ligand_min, dock_obj, dockMethod)
        ligand_min.SetTitle(name)
    else:
        return None, None, None

    res = "{},{},{},{}\n".format(str(pos), name, target_name, score_min)
    return score_min, res, ligand_min

def RunDocking_(smiles: str, dock_obj: oedocking.OEDock, pos: int = None, name: str = None, target_name: str = None, force_flipper: bool = True) -> object:
    """

    :param smiles: a str representing the molecule (SMILES)
    :param dock_obj: The OEDock prepared receptor
    :param pos: int for print string
    :param name: name of ligand for printing string
    :param target_name: name of target from docking
    :param force_flipper: whether or not to flip entamiers when generating conformers
    :return: score, string result, ligand
    """
    if not dock_obj.IsInitialized():
        assert(False)

    scores = []
    ligands = []

    confs = conf_gen.FromString(smiles, force_flipper=force_flipper)
    if confs is None or len(confs) == 0:
        return None, None, None

    for conf in confs:  # dock each Enantiomer
        lig = oechem.OEMol()
        dock_conf.DockConf_(dock_obj, conf, lig, MAX_POSES=1)
        scores.append(dock_obj.ScoreLigand(lig))
        ligands.append(lig)

    # get best score from the Enantiomers
    if len(scores) > 0:
        bs = np.argmin(scores)
        score_min = scores[bs]
        ligand_min = ligands[bs]

        dockMethod = dock_obj.GetName()
        oedocking.OESetSDScore(ligand_min, dock_obj, dockMethod)
        # ligand_min.SetTitle(str(name))
    else:
        return None, None, None

    res = "{},{},{},{},{}\n".format(str(pos), name, smiles, target_name, score_min)
    return score_min, res, ligand_min

def RunDocking_A(smiles, inpath, outpath, target_name,  dock_obj=None, pos=0,
                  name='UNK', oe=False, force_flipper=True, receptor=None):
    from . import conf_gen
    from . import dock_conf
    if not os.path.exists(outpath):
        os.mkdir(outpath)

    confs = conf_gen.FromString(smiles, force_flipper=force_flipper)
    if confs is None or len(confs) == 0:
        return None, None

    scores = []
    ligands = []
    for conf in confs:  # dock each Enantiomer
        lig = oechem.OEMol()
        dock_conf.DockConf_(dock_obj, conf, lig, MAX_POSES=1)
        scores.append(dock_obj.ScoreLigand(lig))
        ligands.append(lig)

    # get best score from the Enantiomers
    if len(scores) > 0:
        bs = np.argmin(scores)
        score_min = scores[bs]
        ligand_min = ligands[bs]

        dockMethod = dock_obj.GetName()
        oedocking.OESetSDScore(ligand_min, dock_obj, dockMethod)
        ligand_min.SetTitle(name)
    else:
        return None, None

    ligand_min.SetTitle(name)
    res = "{},{},{},{},{}\n".format(str(pos), name, smiles, score_min,target_name)


    dock_conf.WriteStructures(receptor, ligand_min, f'{outpath}/apo.pdb', f'{outpath}/lig.pdb', f'{outpath}/com.pdb',
                              oe=oe)
    with open(f'{outpath}/metrics.csv', 'w+') as metrics:
        metrics.write("pos,name,smiles,Dock,receptor\n")
        metrics.write(res)

    return score_min, res


def ParameterizeOE(path):
    """
    Reads in the PDB from 'RunDocking' and outputs 'charged.mol2' of the ligand
    Then runs antechamber to convert this to coordinate (.inpcrd) and 
    parameter (.prmtop) files.
    """
    from openeye import oechem, oequacpac
    mol = oechem.OEMol()
    ifs = oechem.oemolistream()
    if ifs.open(f'{path}/lig.pdb'):
        oechem.OEReadMolecule(ifs, mol)
        ifs.close()
    if not oequacpac.OEAssignCharges(mol, oequacpac.OEAM1BCCCharges()):
        raise (RuntimeError("OEAssignCharges failed."))
    ofs = oechem.oemolostream()
    if ofs.open(f'{path}/charged.mol2'):
        oechem.OEWriteMolecule(ofs, mol)

    import subprocess
    with working_directory(path):
        subprocess.check_output(
            f'antechamber -i lig.pdb -fi pdb -o lig.mol2 -fo mol2 -pf y -an y -a charged.mol2 -fa mol2 -ao crg',
            shell=True)
        subprocess.check_output(f'parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod', shell=True)
        # Wrap tleap
        with open(f'leap.in', 'w+') as leap:
            leap.write("source leaprc.protein.ff14SBonlysc\n")
            leap.write("source leaprc.gaff\n")
            leap.write("set default PBRadii mbondi3\n")
            leap.write("rec = loadPDB apo.pdb # May need full filepath?\n")
            leap.write("saveAmberParm rec apo.prmtop apo.inpcrd\n")
            leap.write("lig = loadmol2 lig.mol2\n")
            leap.write("loadAmberParams lig.frcmod\n")
            leap.write("com = combine {rec lig}\n")
            leap.write("saveAmberParm lig lig.prmtop lig.inpcrd\n")
            leap.write("saveAmberParm com com.prmtop com.inpcrd\n")
            leap.write("quit\n")
        subprocess.check_output(f'tleap -f leap.in', shell=True)


def ParameterizeAMBER(path):
    """
    Alternative method for parameterizing the system. It uses sqm and runs much slower than
    ParameterizeOE, which works through OpenEye.
    This function is pretty much a wrapper for antechamber & tleap. 
    I've kept it as a backup in case there are issues with OpenEye. 
    """
    import subprocess
    with working_directory(path):
        subprocess.check_output(f'antechamber -i lig.pdb -fi pdb -o lig.mol2 -fo mol2 -c bcc -pf y -an y', shell=True)
        subprocess.check_output(f'parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod', shell=True)
        with open(f'leap.in', 'w+') as leap:
            leap.write("source leaprc.protein.ff14SBonlysc\n")
            leap.write("source leaprc.gaff\n")
            leap.write("set default PBRadii mbondi3\n")
            leap.write("rec = loadPDB apo.pdb # May need full filepath?\n")
            leap.write("saveAmberParm rec apo.prmtop apo.inpcrd\n")
            leap.write("lig = loadmol2 lig.mol2\n")
            leap.write("loadAmberParams lig.frcmod\n")
            leap.write("com = combine {rec lig}\n")
            leap.write("saveAmberParm lig lig.prmtop lig.inpcrd\n")
            leap.write("saveAmberParm com com.prmtop com.inpcrd\n")
            leap.write("quit\n")
        subprocess.check_output(f'tleap -f leap.in', shell=True)


def RunMinimization(build_path, outpath, one_traj=False):
    """
    We are minimizing all three structures, then checking the potential energy using GB forcefields
    We could, alternatively, minimize the docked structure and then extract trajectories (1 frame long),
    more like a 1-trajectory mmgbsa.
    output is path used in "RunDocking". It has a metric.csv file.
    """
    from . import minimize
    success = True
    try:
        rec_energy = minimize.MinimizedEnergy(f'{build_path}/apo')
        lig_energy = minimize.MinimizedEnergy(f'{build_path}/lig')
        com_energy = minimize.MinimizedEnergy(f'{build_path}/com')
        diff_energy = com_energy - lig_energy - rec_energy
    except:
        success = False

    if one_traj:
        print("1-traj calculation not ready")
    # TODO: We could decide to do 1-trajectory mmgbsa. It would run about twice as fast as the
    #       current method. I think it would be less accurate, but maybe not. Look into the 1-traj
    #       method from Coveney papers if you want to implement this.

    with open(f'{outpath}/metrics.csv', 'r') as metrics:
        dat = metrics.readlines()
    with open(f'{outpath}/metrics.csv', 'w') as metrics:
        metrics.write(dat[0].replace('\n', ',Minimize,Minimize_U\n'))
        if success:
            metrics.write(dat[1].replace('\n', ',{},{}\n'.format(diff_energy, 0)))
        else:
            metrics.write(dat[1].replace('\n', ',NA,NA\n'))


def RunMinimization_(build_path, outpath, one_traj=False, write=False, gpu=False):
    from . import minimize
    success = True
    try:
        lig_energy = minimize.MinimizedEnergy(f'{build_path}/lig', gpu=gpu)
        print('lig', lig_energy)
        rec_energy = minimize.MinimizedEnergy(f'{build_path}/apo', gpu=gpu)
        print('rec_energy', rec_energy)
        com_energy = minimize.MinimizedEnergy(f'{build_path}/com', gpu=gpu)
        print('com_energy', com_energy)

        diff_energy = com_energy - lig_energy - rec_energy
    except:
        success = False

    if one_traj:
        print("1-traj calculation not ready")
    # TODO: We could decide to do 1-trajectory mmgbsa. It would run about twice as fast as the
    #       current method. I think it would be less accurate, but maybe not. Look into the 1-traj
    #       method from Coveney papers if you want to implement this.
    if write:
        with open(f'{outpath}/metrics.csv', 'r') as metrics:
            dat = metrics.readlines()
        with open(f'{outpath}/metrics.csv', 'w') as metrics:
            metrics.write(dat[0].replace('\n', ',Minimize,Minimize_U\n'))
            if success:
                metrics.write(dat[1].replace('\n', ',{},{}\n'.format(diff_energy, 0)))
            else:
                metrics.write(dat[1].replace('\n', ',NA,NA\n'))
    if success:
        return diff_energy
    else:
        return np.nan


def Simulation_explicit(inpath, outpath, nsteps, comp='com'):
    if not os.path.exists(outpath):
        os.mkdir(outpath)

    from . import minimize
    success = True
    try:
        potential = minimize.simulation(f'{inpath}/{comp}', outpath, nsteps)
    except:
        success = False

    with open(f'{inpath}/metrics.csv', 'r') as metrics:
        dat = metrics.readlines()
    with open(f'{inpath}/metrics.csv', 'w') as metrics:
        metrics.write(dat[0].replace('\n', ',U_minimized_explicit\n'))
        if success:
            metrics.write(dat[1].replace('\n', ',{}\n'.format(potential)))
        else:
            metrics.write(dat[1].replace('\n', ',NA\n'))
    if success:
        return potential
    else:
        return np.nan


def RunMMGBSA(inpath, outpath, niter=1000):
    """
    1 'iteration' corresponds to 1 ps.
    """
    from . import mmgbsa
    crds = {'lig': f'{inpath}/lig.inpcrd', 'apo': f'{inpath}/apo.inpcrd', 'com': f'{inpath}/com.inpcrd'}
    prms = {'lig': f'{inpath}/lig.prmtop', 'apo': f'{inpath}/apo.prmtop', 'com': f'{inpath}/com.prmtop'}

    enthalpies = mmgbsa.simulate(crds, prms, niter)
    # enthalpies is a list of energies from each iteration
    mmgbsa.subsample(enthalpies)
    # We subsample the enthalpies using a method from John Chodera that determines the equilibration
    #   and autocorrelation times. This allows us to extract an uncertainty.
    #   See the file mmgbsa.py or his package 'pymbar' for more detail.
    energies = mmgbsa.mmgbsa(enthalpies)

    with open(f'{outpath}/metrics.csv', 'r') as metrics:
        dat = metrics.readlines()
    with open(f'{outpath}/metrics.csv', 'w') as metrics:
        metrics.write(dat[0].replace('\n', ',mmgbsa,mmgbsa_U\n'))
        metrics.write(dat[1].replace('\n', ',{},{}\n'.format(energies[0]['diff'], energies[1]['diff'])))
    return energies


def RunMMGBSA_(inpath, outpath, gpu=False, niter=1000):
    """
    1 'iteration' corresponds to 1 ps.
    """
    from . import mmgbsa
    crds = {'lig': f'{inpath}/lig.inpcrd', 'apo': f'{inpath}/apo.inpcrd', 'com': f'{inpath}/com.inpcrd'}
    prms = {'lig': f'{inpath}/lig.prmtop', 'apo': f'{inpath}/apo.prmtop', 'com': f'{inpath}/com.prmtop'}

    enthalpies = mmgbsa.simulate(crds, prms, outpath, gpu=gpu, niters=niter)
    # enthalpies is a list of energies from each iteration
    mmgbsa.subsample(enthalpies)
    # We subsample the enthalpies using a method from John Chodera that determines the equilibration
    #   and autocorrelation times. This allows us to extract an uncertainty.
    #   See the file mmgbsa.py or his package 'pymbar' for more detail.
    energies = mmgbsa.mmgbsa(enthalpies)

    with open(f'{outpath}/metrics.csv', 'r') as metrics:
        dat = metrics.readlines()
    with open(f'{outpath}/metrics.csv', 'w') as metrics:
        metrics.write(dat[0].replace('\n', ',mmgbsa,mmgbsa_U\n'))
        metrics.write(dat[1].replace('\n', ',{},{}\n'.format(energies[0]['diff'], energies[1]['diff'])))
    return energies[0]['diff']


def RunAlchemy(path, niter=2500, nsteps_per_iter=1000, nlambda=11):
    """
    Default is a 5 ns simulation with sampling every 2 ps
    """
    from . import alchemy
    [energy, err] = alchemy.SimulateAlchemy(path, niter, nsteps_per_iter, nlambda)
    with open(f'{path}/metrics.csv', 'r') as metrics:
        dat = metrics.readlines()
    with open(f'{path}/metrics.csv', 'w') as metrics:
        metrics.write(dat[0].replace('\n', ',alchemy,alchemy_U\n'))
        metrics.write(dat[1].replace('\n', ',{},{}\n'.format(energy, err)))
    return energy, err
