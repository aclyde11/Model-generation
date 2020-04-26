

# Docking
There are four flavors of docking scripts here.
- docking.py is a script aimed at running on a small job on a local laptop or computer. It can be run in parallel on a single machine using python multiprocressing. Not suitable for HPC.
- theta_dock.py is a script aimed at HPC. It does not use any load balancing and simply splits the input file over the ranks. Does not need MPI4py.
```shell script
aprun -n <number of ranks> ... python theta_dock.py -i <input.csv> -o <out sdf locations: outfolder/out.sdf> -v -n <number of ranks> -r <receptor file>
cat outfolder/*.sdf > all.sdf
python convert.py -i all.sdf -o all_cleaned.sdf -u
```
- theta_dock_mpi.py is a script aimed at HPC. It uses mpi4py to do load balancing while the job is running, handing out data in chunks of 5.
```shell script
aprun -n <number of ranks> ... python theta_dock.py -i <input.csv> -o <out sdf locations: outfolder/out.sdf> -v -r <receptor file>
cat outfolder/*.sdf > all.sdf
python convert.py -i all.sdf -o all_cleaned.sdf -u
```
- theta_dock_mpi_conf.py is a script aimed at HPC. it uses mpi4py to do load balancing. It also uses a conformer oeb input file to skip generating conformers every time resulting in significant speed up. The data requires some post processing however. 
```shell script
aprun -n <number of ranks> ... python theta_dock_mpi_conf.py -i <input.oeb> -o <out sdf locations: outfolder/out.sdf> -v -r <receptor file>
cat outfolder/*.sdf > all.sdf
python convert.py -i all.sdf -o all_cleaned.sdf -u
```

While these scripts are running, one can just query the output to see the progress (i.e. wc -l 4303.output). Except the theta_dock_mpi_conf will not show correctly, one needs to use python checkcount.py <outputlog.output> to get the correct count of moleculees procressed so far. 


## Helper Scripts
There are three utility scripts ```convert.py, sdsort.py, enamine_fix_oe_title.py```
All scripts use openeye for IO, so input and output formats are flexible and open-ended (-i/o {.sdf, .smi, .csv, .mol2, .pdb}). Openeye sets all titles correctly and calls the names in CSV TITLE and SMILES. Docking scores can appear as "[FRED, HYBRID, ''] Chemguass4." Autodock rescores get tagged AutodockVinaRescoreOnly. There is an issue where 

- ```sdsort.py```
This script performs a sort based an SD tag. For example, to sort an SDF file by Chemgauss one can do
```shell script
python sdsort.py -i out.sdf -o out_sorted.sdf -t "Chemgauss4" -u
```
where one can add the ```-u``` to make sure some molecules are aggregated by title as Omega during the conformer generation may generate a different molecular entry, even though it originated from the same SMILES.

You can convert in one step and optionally only get the top n as well:
```shell script
python sdsort.py -i out.sdf -o out_sorted_top100.csv -t "Chemgauss4" -u -n 100
```

- ```convert.py```
This script performs basic conversion.
```shell script
python convert.py -i out.sdf -o out_sorted.csv
```
where one can add the ```-u``` to make sure some molecules are aggregated by title as Omega during the conformer generation may generate a different molecular entry, even though it originated from the same SMILES.

- ```enamine_fix_oe_title.py```
This script fixes a ~bug~ where Enamine molecules sometimes don't get titled but rather their name goes to an SD tag "Catalog ID." This finds those molecules and fixes the bug.