from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit, openmm
from openforcefield.topology import Molecule
from openmmforcefields.generators import SystemGenerator


def MinimizedEnergyGAFF(filepath, ligand_path, cache, gpu=False):
    water_model = 'tip3p'
    solvent_padding = 10.0 * unit.angstrom
    box_size = openmm.vec3.Vec3(3.4, 3.4, 3.4) * unit.nanometers
    ionic_strength = 100 * unit.millimolar  # 100
    pressure = 1.0 * unit.atmospheres
    collision_rate = 91.0 / unit.picoseconds
    temperature = 310.15 * unit.kelvin
    timestep = 4.0 * unit.femtoseconds
    nsteps_per_iteration = 250
    iterations = 1000
    protein_forcefield = 'amber14/protein.ff14SB.xml'
    small_molecule_forcefield = 'gaff-2.11' # only if you really like atomtypes
    solvation_forcefield = 'amber14/tip3p.xml'


    off_molecule = Molecule.from_file(f'{ligand_path}.pdb')

    barostat = openmm.MonteCarloBarostat(pressure, temperature)

    common_kwargs = {'removeCMMotion': True, 'ewaldErrorTolerance': 5e-04,
                     'nonbondedMethod': app.PME, 'hydrogenMass': 3.0 * unit.amu}
    # unconstrained_kwargs = {'constraints': None, 'rigidWater': False}
    constrained_kwargs = {'constraints': app.HBonds, 'rigidWater': True}
    forcefields = [protein_forcefield, solvation_forcefield]

    openmm_system_generator = SystemGenerator(forcefields=forcefields,
                                              molecules=[off_molecule],
                                              small_molecule_forcefield=small_molecule_forcefield, cache=cache,
                                              barostat=barostat,
                                              forcefield_kwargs={**common_kwargs, **constrained_kwargs})

    with open(f'{filepath}.pdb', 'r') as infile:
        lines = [line for line in infile if 'UNK' not in line]
    from io import StringIO
    pdbfile_stringio = StringIO(''.join(lines))

    # Read the unsolvated system into an OpenMM Topology
    pdbfile = app.PDBFile(pdbfile_stringio)
    topology, positions = pdbfile.topology, pdbfile.positions

    print('Adding solvent...')
    modeller = app.Modeller(topology, positions)
    kwargs = {'padding': solvent_padding}
    modeller.addHydrogens(openmm_system_generator.forcefield)
    modeller.addSolvent(openmm_system_generator.forcefield, model='tip3p', ionicStrength=ionic_strength, **kwargs)
    print('Building system...')

    # Create an OpenMM system
    system = openmm_system_generator.create_system(modeller.topology)
    print("system made")
    if gpu:
        platform = openmm.Platform.getPlatformByName('CUDA')
        platform.setPropertyDefaultValue('Precision', 'mixed')
    else:
        platform = openmm.Platform.getPlatformByName('CPU')

    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
    integrator.setConstraintTolerance(0.00001)
    simulation = app.Simulation(modeller.topology, system, integrator, platform=platform)
    simulation.context.setPositions(modeller.positions)
    print("wr")

    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(310.15 * unit.kelvin)
    simulation.step(100)
    simulation.reporters.append(app.PDBReporter(f'{filepath}_output.pdb', 100))

    simulation.step(5000)

    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule / unit.mole)
    with open(f'{filepath}.state', 'w') as f:
        simulation.saveState(f)
    return energy


def MinimizedEnergy(filepath, gpu=False):
    prmtop = app.AmberPrmtopFile(f'{filepath}.prmtop')
    inpcrd = app.AmberInpcrdFile(f'{filepath}.inpcrd')
    print(f'{filepath}.prmtop')
    system = prmtop.createSystem(implicitSolvent=app.GBn2,
                                 nonbondedMethod=app.CutoffNonPeriodic,
                                 nonbondedCutoff=1.0*unit.nanometers,
                                 constraints=app.HBonds,
                                 rigidWater=True,
                                 ewaldErrorTolerance=0.0005)

    integrator = mm.LangevinIntegrator(310.15*unit.kelvin,
                                       1.0/unit.picoseconds,
                                       2.0*unit.femtoseconds)
    integrator.setConstraintTolerance(0.00001)
    # TODO: This should just recognize whatever the computer is capable of, not force CUDA.

    if gpu:
        platform = 'CUDA'
    else:
        platform = 'CPU'
    platform = mm.Platform.getPlatformByName(platform)

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

def simulation(filepath, outpath, nsteps, gpu=True):
    prmtop = app.AmberPrmtopFile(f'{filepath}.prmtop')
    inpcrd = app.AmberInpcrdFile(f'{filepath}.inpcrd')
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    modeller = app.Modeller(prmtop.topology, inpcrd.positions)
    modeller.addSolvent(forcefield, padding=1.4*unit.nanometer)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometer,
            constraints=app.HBonds)
    integrator = mm.LangevinIntegrator(310.15*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picosecond)
    if gpu:
        platform = 'CUDA'
    else:
        platform = 'CPU'
    platform = mm.Platform.getPlatformByName(platform)
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
    
