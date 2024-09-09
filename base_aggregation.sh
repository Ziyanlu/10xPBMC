#!/bin/bash
#SBATCH -A b1042 ## EDIT THIS TO BE YOUR ALLOCATION 
#SBATCH -p genomics ## <-- EDIT THIS TO BE YOUR QUEUE NAME
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=48G
#SBATCH --job-name=10xGenomic_data
#SBATCH --output=outlog
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ziyan@northwestern.edu

module purge all                ## Unload existing modules
source /projects/b1038/Pulmonary/sfenske/projects/reference_atlas/code/ref_atlas_venv/bin/activate # activate virtual environment with 
# packages I need (scvi, scanpy)

deactivate

#python /projects/b1038/Pulmonary/ziyan/10xGenomics_PBMC_datasets/code/base_aggregation.py     ## Run the program
#sample_adata = sc.read_10x_h5(/projects/b1038/Pulmonary/ziyan/10xGenomics_PBMC_datasets/raw/5p_nextgem_5K_PBMC/filtered_feature_bc_matrix.h5)