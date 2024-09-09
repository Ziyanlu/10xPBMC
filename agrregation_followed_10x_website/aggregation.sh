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
module load bcl2fastq/2.20      ## load necessary modules (software, libraries)

/projects/b1038/tools/cellranger-7.2.0/cellranger aggr \
--id=agg10 \
--csv=aggregation_PBMC.csv \
--output-dir=/projects/b1038/Pulmonary/ziyan/10xGenomics_PBMC_datasets/data
