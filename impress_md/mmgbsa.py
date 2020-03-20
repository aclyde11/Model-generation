import numpy as np
from pymbar import timeseries
from math import sqrt


def simulate(inpcrd_filenames, prmtop_filenames, path, niterations=5000, implicit=True, gpu=False):
    """
    The program simulates three systems: the ligand alone, protein alone, and complex.
    Input is a dict of files to the input coordinates (.inpcrd) and parameters (.prmtop) 
      as well as the number of iterations. One iteration is one picosecond. 
    Output is a dict of a list of the ennthalpies calculated using mmgbsa for each system.
    """
    from simtk.openmm import app
    import simtk.openmm as mm
    from simtk import unit
    
    phases = inpcrd_filenames.keys()
    nsteps_per_iteration = 500 # 1 picosecond
    
    enthalpies = dict()
    for phase in phases:
        print("on phase", phase)
        enthalpies[phase] = np.zeros([niterations])
        prmtop = app.AmberPrmtopFile(prmtop_filenames[phase])
        inpcrd = app.AmberInpcrdFile(inpcrd_filenames[phase])
        
        system = prmtop.createSystem(implicitSolvent=app.GBn2, 
                 nonbondedMethod=app.CutoffNonPeriodic,
                 nonbondedCutoff=2.0*unit.nanometers, 
                 constraints=app.HBonds)
        integrator = mm.LangevinIntegrator(310.15*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
        integrator.setConstraintTolerance(0.00001)
        
        platform = mm.Platform.getPlatformByName('CUDA' if gpu else 'CPU')
        # I'm not sure what the CudaDeviceIndex is about. Austin at 
        properties = {'CudaPrecision': 'mixed', 'CudaDeviceIndex' : '0'}
        simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)
        simulation.context.setPositions(inpcrd.positions)
        if phase == 'com':
            simulation.reporters.append(app.PDBReporter(path + '/' + phase + '_output.dcd', 10000))
        # Minimize & equilibrate
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(310.15*unit.kelvin)
        simulation.step(100)

        # Run simulation
        for iteration in range(niterations):
            simulation.step(nsteps_per_iteration)
            state = simulation.context.getState(getEnergy=True)
            potential_energy = state.getPotentialEnergy()
            enthalpies[phase][iteration] = potential_energy.value_in_unit(unit.kilojoules_per_mole)
        del simulation
        del system
        del platform
    return enthalpies


def subsample(enthalpies):
    """
    Subsamples the enthalpies using John Chodera's code.
    This is probably better than the simple cutoff we normally use.
    No output -- it modifies the lists directly
    """
    for phase in enthalpies:
        # Use automatic equilibration detection and pymbar.timeseries to subsample
        [t0, g, Neff_max] = timeseries.detectEquilibration(enthalpies[phase])
        enthalpies[phase] = enthalpies[phase][t0:]
        indices = timeseries.subsampleCorrelatedData(enthalpies[phase], g=g)
        enthalpies[phase] = enthalpies[phase][indices]


def mmgbsa(enthalpies):
    """
    Returns DeltaG, errDeltaG : float
        Estimated free energy of binding
    """
    DeltaH = dict()
    varDeltaH = dict()
    errDeltaH = dict()
    # TODO: it calculates standard error rather than variance
    for phase in enthalpies:
        DeltaH[phase] = enthalpies[phase].mean()
        varDeltaH[phase] = enthalpies[phase].std()**2
        errDeltaH[phase] = varDeltaH[phase]/len(enthalpies[phase])

    # TODO: Shouldn't have to use names?
    try:
        DeltaH['diff'] = 2*DeltaH['complex']
    except:
        DeltaH['diff'] = 2*DeltaH['com']
    varDeltaH['diff'] = 0
    errDeltaH['diff'] = 0
    for phase in enthalpies:
        DeltaH['diff'] -= DeltaH[phase]
        varDeltaH['diff'] += varDeltaH[phase]
        errDeltaH['diff'] += errDeltaH[phase]

    errDeltaH['diff'] = sqrt(errDeltaH['diff'])
    return DeltaH, errDeltaH

