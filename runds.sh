#!/bin/bash
#COBALT -t 30
#COBALT -n 4
#COBALT -q debug-cache-quad
#COBALT -A CVD_Research

module load  miniconda-3.6/conda-4.5.12
source activate oedock
export OE_LICENSE=/home/aclyde/oe_license.txt 
export PATH=/projects/candle_aesp/aclyde/amber18/bin:$PATH
source /projects/candle_aesp/aclyde/amber18/amber.sh

cd /projects/candle_aesp/aclyde/Model-generation
aprun -n 1024 -N 256 -j 4 -d 1  python runner.py --smiles input/enamine_diverse.smi --receptor_file input/swiss_plpro.oeb --path swis_plpro_test --dock_only --target_name plpro_swiss --dbase_name test
