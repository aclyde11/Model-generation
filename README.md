

# Docking
 The smiles input file should have a column that matches something like SMILES or smiles, and a name column that matches something like ID, CID, TITLE, NAME etc (case doesn't matter). -n will control number of workers to run, default is 1 (single threaded). -v will print somethings to screen while it runs. The output file format is automatically detected. SDF recomended. Use convert.py or sdsorter.py to sort or manage that conversion after docking.
```shell script
python docking.py -i my_smiles.csv -o results.{sdf, csv} -n 1 -r recetor.oeb -v 
```

```shell script
python convert.py -i input.sdf -o output.csv
```

```shell script
python sdsort.py -i input.sdf -o output.sdf -t "FRED Chemgauss4 Score"
```

#### Don't rec venturing below this

# COVID WORK TO DO 

Hello,

Major work items are creating a robust workflow for the giga-docking. Please make sure you're on the vcovid branch


Look at the function in the inferface_function file called Run_Docking_
The arguments are as follows:

There are two options to run.

Preload the docking receptor 

smiles -- a string of the smiles
inpath -- a receptor file which will be provided as oeb
dbase_name -- a string which is the name of the dbase run where the smiles comes
target_name -- a string which is the name of the target 
pos=0 --file string index position of the thing
receptor_file=None -- a receptor file which will be provided as oeb 
dock_obj=None -- the dock obj is you are precomputing it
recept=None -- the receptor object if you are precomputing it
name='UNK'
docking_only=False
 
Do not preload the docking receptor 
smiles -- a string of the smiles
inpath -- a receptor file which will be provided as oeb
dbase_name -- a string which is the name of the dbase run where the smiles comes
target_name -- a string which is the name of the target 
pos=0 -- file string index position of the thing
receptor_file=None -- a receptor file which will be provided as oeb 
dock_obj=None -- the dock obj is you are precomputing it
recept=None -- the receptor object if you are precomputing it
name='UNK'
docking_only=False


## INTERFACE

a list of receptor file .oeb, and a csv of smiles. 

receptor file on one axis
smiles on the other. 


So in the folder input, take the enamine_diverse.smi and the swiss_plpro.oeb. Test on those. Sorry if my recent changes broke this code, but should be easy to fix.

Taake a look at theta dock to see how this works.
---------

This runs of x86 ONLY for the moment. 

Build the conda enviroment using the yml file or just run stuff and figure out how to make it work (it's not rocket science ;)

You can run end to end with 
```python 
mpiexec -np 8 python runner.py input/john_smiles_kinasei.smi
```

You can alter how those runs peform using the policy file.

In the future, this workflow will be moved to Radical/ENTK, but for now this works. 



---- 

# Model-generation
Python scripts to generate an MD-ready model from smiles strings and run simple free energy calculations 

There are two command line executables.

## docking.py
* Takes a smiles and pdb, generates conformers, docks, and scores the ligand.
* The output is a set of simulation-ready structures (ligand, apo, and complex) and a file called metrics.csv, which has the docking score and associated uncertainties. Most uncertainties are 0 right now. There are other auxiliary files that are saved in the output directory.
* Dependencies: OpenEye, Ambertools, Ambermini, docopt

## param.py
* Parameterizes the ligand using either OpenEye or Amber.

## mmgbsa.py
* This command should only be run after docking.py. 
* Input is a path to the directory where the input coordinates and parameters are saved. This should be the output path from the docking.py command.
* Also takes the nanosecond length of the simulation. 0 corresponds to an energy minimization.
* Output adds to the metrics.csv file
* Dependencies: OpenMM, numpy, pymbar, docopt

## alchem.py
* Uses an alchemical method to calculate the absolute binding free energy of a ligand.
* User enters the number of lambda windows and the length of simulation at each window.
* The windows run in series but this could be parallelized.
* Applying constraints should increase the convergence of the system

To get the four metrics for a smiles, including a 5 ns simulation, pick a smiles and call
~~~bash
python docking.py -s $SMILES -i "input/receptor.oeb" -o "test"
python param.py -i "test"
python mmgbsa.py -p "test" -n 0
python mmgbsa.py -p "test" -n 5
python alchem.py -i "test" -l 6 -n 5
~~~
