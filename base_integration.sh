#!/bin/bash
#SBATCH -A b1042 ## EDIT THIS TO BE YOUR ALLOCATION 
#SBATCH -p genomics ## <-- EDIT THIS TO BE YOUR QUEUE NAME
#SBATCH --time=12:00:00
#SBATCH --job-name=10xGenomic_data
#SBATCH --output=outlog.txt
#SBATCH --error=errlog
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ziyan@northwestern.edu

module purge all                                                                                   ## Unload existing modules

source /projects/b1038/Pulmonary/sfenske/projects/reference_atlas/code/ref_atlas_venv/bin/activate # activate virtual environment with 
# packages I need (scvi, scanpy)
python /projects/b1038/Pulmonary/ziyan/10xGenomics_PBMC_datasets/code/v1_integration_from_Sam.py   #run the code

deactivate

