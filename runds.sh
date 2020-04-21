#!/bin/bash
#COBALT -t 12:00
#COBALT -n 64
#COBALT -q R.CVD_Research
#COBALT -A CVD_Research

module load  miniconda-3.6/conda-4.5.12
source activate oedock
export OE_LICENSE=/home/aclyde/oe_license.txt 
export PATH=/projects/candle_aesp/aclyde/amber18/bin:$PATH
source /projects/candle_aesp/aclyde/amber18/amber.sh

cd /projects/candle_aesp/aclyde/Model-generation


unset PYTHONPATH
unset LD_LIBRARY_PATH

module load datascience/mpi4py
export PYTHONPATH=$PYTHONPATH:/gpfs/mira-home/aclyde/.conda/envs/oedock/lib/python3.6/site-packages/
export LD_LIBRARY_PATH=/gpfs/mira-home/aclyde/.conda/envs/oedock/lib/:$LD_LIBRARY_PATH
export ATP_ENABLED=1

aprun -n 8192 -N 128 -j 2 -d 1  python runner.py --smiles input/pubchem_compounds.smi --receptor_file input/swiss_plpro.oeb --path swis_plpro_pcc --dock_only --target_name plpro_swiss --dbase_name pubchemcompounds
