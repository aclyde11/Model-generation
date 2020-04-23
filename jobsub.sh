#!/bin/bash
#COBALT -t 4:00 
#COBALT  -n 96
#COBALT -q CVD_Research
#COBALT -A CVD_Research

echo "Starting Cobalt job script on 128 nodes with 64 ranks on each node" 
echo "with 4 OpenMP threads per rank, for a total of 8192 ranks and 32768 OpenMP threads" 
module load miniconda-3
conda activate /lus/theta-fs0/projects/CVD_Research/aclyde/openeyeenv
cd /projects/CVD_Research/aclyde/run1
export OE_LICENSE=/home/aclyde/oe_license.txt
aprun -n 24576 -N 256 -cc depth -d 1 -j 4 python ../Model-generation/theta_dock.py -i /projects/CVD_Research/aclyde/Workflow0COVID/smiles/BL2.csv  -o out.sdf -v -n 24576 -r /projects/CVD_Research/aclyde/Workflow0COVID/receptors/receptorsVA/ADRP_ADPR_pocket1_rec_6W02_A__DU__APR_A-201.oeb.gz
