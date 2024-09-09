#!/bin/bash
#SBATCH -A b1042 ## EDIT THIS TO BE YOUR ALLOCATION 
#SBATCH -p genomics ## <-- EDIT THIS TO BE YOUR QUEUE NAME
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=48G
#SBATCH --job-name=10xGenomic_data
#SBATCH --output=outlog
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ziyan@northwestern.edu

module purge all                ## Unload existing modules
module load bcl2fastq/2.20      ## load necessary modules (software, libraries)

# Output Files
curl -O https://cf.10xgenomics.com/samples/cell-exp/8.0.0/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx_web_summary.html
curl -O https://cf.10xgenomics.com/samples/cell-exp/8.0.0/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx_metrics_summary.csv
curl -O https://cf.10xgenomics.com/samples/cell-exp/8.0.0/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx_count_sample_cloupe.cloupe
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/8.0.0/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx_count_sample_alignments.bam
curl -O https://cf.10xgenomics.com/samples/cell-exp/8.0.0/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx_count_sample_alignments.bam.bai
curl -O https://cf.10xgenomics.com/samples/cell-exp/8.0.0/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx_count_sample_filtered_barcodes.csv
curl -O https://cf.10xgenomics.com/samples/cell-exp/8.0.0/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx_count_sample_filtered_feature_bc_matrix.h5
curl -O https://cf.10xgenomics.com/samples/cell-exp/8.0.0/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx_count_sample_filtered_feature_bc_matrix.tar.gz
curl -O https://cf.10xgenomics.com/samples/cell-exp/8.0.0/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx_count_sample_molecule_info.h5
curl -O https://cf.10xgenomics.com/samples/cell-exp/8.0.0/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx_count_feature_reference.csv
curl -O https://cf.10xgenomics.com/samples/cell-exp/8.0.0/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx_count_analysis.tar.gz
