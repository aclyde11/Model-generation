import sys, os
from contextlib import contextmanager

@contextmanager
def working_directory(directory):
    owd = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)

def RunDocking(smiles, inpath, outpath, padding=4, return_scores=False, write_metrics_out=True):
    from . import conf_gen
    from . import dock_conf
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    confs = conf_gen.SelectEnantiomer(conf_gen.FromString(smiles))
    # This receptor can be pre-compiled to an oeb. It may speed things up notably
    filename, file_extension = os.path.splitext(inpath)
    if file_extension == ".oeb":
        receptor = dock_conf.PrepareReceptorFromBinary(inpath)
    else: # else it is a pdb
        receptor = dock_conf.PrepareReceptor(inpath,padding,outpath)
    dock, lig = dock_conf.DockConf(receptor,confs,MAX_POSES=1)

    dock_conf.WriteStructures(receptor, lig, f'{outpath}/apo.pdb', f'{outpath}/lig.pdb')
    best_d_score = dock_conf.BestDockScore(dock,lig)
    if write_metrics_out:
        with open(f'{outpath}/metrics.csv','w+') as metrics:
            metrics.write("Dock,Dock_U\n")
            metrics.write("{},{}\n".format(best_d_score,0))
    # from openeye import oedepict
    # oedepict.OEPrepareDepiction(lig)
    # oedepict.OERenderMolecule(f'{outpath}/lig.png',lig)
    if return_scores:
        return best_d_score

def ParameterizeOE(path):
    """
    Reads in the PDB from 'RunDocking' and outputs 'charged.mol2' of the ligand
    """
    from openeye import oechem, oeomega, oequacpac
    mol = oechem.OEMol()
    ifs = oechem.oemolistream()
    if ifs.open(f'{path}/lig.pdb'):
        oechem.OEReadMolecule(ifs,mol)
        ifs.close()
    if not oequacpac.OEAssignCharges(mol,oequacpac.OEAM1BCCCharges()):
        print("Error parameteririzng ", path)
        return
    ofs = oechem.oemolostream()
    if ofs.open(f'{path}/charged.mol2'):
        oechem.OEWriteMolecule(ofs,mol)
    
    import subprocess
    with working_directory(path):
        subprocess.check_output(f'antechamber -i lig.pdb -fi pdb -o lig.mol2 -fo mol2 -pf y -an y -a charged.mol2 -fa mol2 -ao crg',shell=False)
        subprocess.check_output(f'parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod',shell=False)
        # Wrap tleap to get $path/(com|lig|apo).(inpcrd|prmtop)
        with open(f'leap.in','w+') as leap:
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
        subprocess.check_output(f'tleap -f leap.in',shell=True)


# DEPRICATED -- this doesn't run on rhea properly for some reason
# It may be due to a write-permission on Rhea that is now fixed.
# In any event, I'm running it in OE now.
def ParameterizeSystem(path):
    import subprocess
    with working_directory(path):
        subprocess.check_output(f'antechamber -i lig.pdb -fi pdb -o lig.mol2 -fo mol2 -c bcc -pf y -an y',shell=True)
        subprocess.check_output(f'parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod',shell=True)
        # Wrap tleap to get $path/(com|lig|apo).(inpcrd|prmtop)
        with open(f'leap.in','w+') as leap:
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
        subprocess.check_output(f'tleap -f leap.in',shell=True)
    

def RunMinimization(build_path, outpath, one_traj=False, return_score=False, write_metrics_out=True):
    """
    We are minimizing all three structures, then checking the potential energy using GB forcefields
    We could, alternatively, minimize the docked structure and then extract trajectories (1 frame long), more like a 1-trajectory mmgbsa.
    output is path used in "RunDocking". It has a metric.csv file.
    Notes
        MinimizedEnergy function requires path/prefix to inpcrd and prmtop files.
        These are created by the ParameterizeSystem function...
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

    if write_metrics_out:
        with open(f'{outpath}/metrics.csv','r') as metrics:
            dat = metrics.readlines()
        with open(f'{outpath}/metrics.csv','w') as metrics:
            metrics.write(dat[0].replace('\n',',Minimize,Minimize_U\n'))
            if success:
                metrics.write(dat[1].replace('\n',',{},{}\n'.format(diff_energy,0)))
            else:
                metrics.write(dat[1].replace('\n',',NA,NA\n'))
    if return_score:
        return diff_energy

def RunMMGBSA(inpath, outpath, niter=1000, write_metrics_out=True, return_score=False):
    from . import mmgbsa_new as mmgbsa
    crds = {'lig':f'{inpath}/lig.inpcrd','apo':f'{inpath}/apo.inpcrd','com':f'{inpath}/com.inpcrd'}
    prms = {'lig':f'{inpath}/lig.prmtop','apo':f'{inpath}/apo.prmtop','com':f'{inpath}/com.prmtop'}
    # print("Starting simulation...")
    enthalpies = mmgbsa.simulate(crds, prms, niter)
    # print("Subsampling to reduce variance")
    mmgbsa.subsample(enthalpies)
    energies = mmgbsa.mmgbsa(enthalpies)
    # Now write this to file
    if write_metrics_out:
        with open(f'{outpath}/metrics.csv','r') as metrics:
            dat = metrics.readlines()
        with open(f'{outpath}/metrics.csv','w') as metrics:
            metrics.write(dat[0].replace('\n',',mmgbsa,mmgbsa_U\n'))
            metrics.write(dat[1].replace('\n',',{},{}\n'.format(energies[0]['diff'],energies[1]['diff'])))

    if return_score:
        return energies[0]['diff'],energies[1]['diff']
    else:
        return energies
