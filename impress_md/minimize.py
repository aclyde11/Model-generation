from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit

def MinimizedEnergy(filepath):
    prmtop = app.AmberPrmtopFile(f'{filepath}.prmtop')
    inpcrd = app.AmberInpcrdFile(f'{filepath}.inpcrd')
    system = prmtop.createSystem(implicitSolvent=app.GBn2,
                                 nonbondedMethod=app.CutoffNonPeriodic,
                                 nonbondedCutoff=1.0*unit.nanometers,
                                 constraints=app.HBonds,
                                 rigidWater=True,
                                 ewaldErrorTolerance=0.0005)

    integrator = mm.LangevinIntegrator(300*unit.kelvin,
                                       1.0/unit.picoseconds,
                                       2.0*unit.femtoseconds)
    integrator.setConstraintTolerance(0.00001)
    # TODO: This should just recognize whatever the computer is capable of, not force CUDA.
    platform = mm.Platform.getPlatformByName('CPU')
    # TODO: I am not sure if mixed precision is necessary. It dramatically changes the results.

    simulation = app.Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(inpcrd.positions)
    
    simulation.minimizeEnergy()
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule/unit.mole)
    return energy


# def MinimizedEnergyWithParam(filepath):
#     prmtop = app.AmberPrmtopFile(f'{filepath}.prmtop')
#     inpcrd = app.AmberInpcrdFile(f'{filepath}.inpcrd')
#     system = prmtop.createSystem(implicitSolvent=app.GBn2,
#                                  nonbondedMethod=app.CutoffNonPeriodic,
#                                  nonbondedCutoff=1.0 * unit.nanometers,
#                                  constraints=app.HBonds,
#                                  rigidWater=True,
#                                  ewaldErrorTolerance=0.0005)
#
#     integrator = mm.LangevinIntegrator(300 * unit.kelvin,
#                                        1.0 / unit.picoseconds,
#                                        2.0 * unit.femtoseconds)
#     integrator.setConstraintTolerance(0.00001)
#     # TODO: This should just recognize whatever the computer is capable of, not force CUDA.
#     platform = mm.Platform.getPlatformByName('CUDA')
#     # TODO: I am not sure if mixed precision is necessary. It dramatically changes the results.
#     properties = {'CudaPrecision': 'mixed'}
#
#     simulation = app.Simulation(prmtop.topology, system, integrator, platform)
#     simulation.context.setPositions(inpcrd.positions)
#
#     simulation.minimizeEnergy()
#     energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule / unit.mole)
#     return energy

def simulation(filepath, outpath, nsteps):
    prmtop = app.AmberPrmtopFile(f'{filepath}.prmtop')
    inpcrd = app.AmberInpcrdFile(f'{filepath}.inpcrd')
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    modeller = app.Modeller(prmtop.topology, inpcrd.positions)
    modeller.addSolvent(forcefield, padding=1.4*unit.nanometer)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometer,
            constraints=app.HBonds)
    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picosecond)
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'double'}
    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    if nsteps != 0:
        simulation.reporters.append(app.DCDReporter(f'{outpath}/traj.dcd', 25000))
        simulation.reporters.append(app.StateDataReporter(f'{outpath}/sim.log', 25000, step=True,
        potentialEnergy=True, temperature=True))
        simulation.reporters.append(app.CheckpointReporter(f'{outpath}/traj.chk', 250000))
        simulation.step(nsteps)
        positions = simulation.context.getState(getPositions=True).getPositions()
        app.PDBFile.writeFile(simulation.topology, positions, open(f'{outpath}/output.pdb', 'w'))
    potential = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule/unit.mole)
    return potential
    
