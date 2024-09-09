import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import warnings
import scanpy as sc
import scvi
import sc_utils
import scrublet as scr
warnings.filterwarnings('ignore')
    
cr_dir = '/projects/b1038/Pulmonary/ziyan/10xGenomics_PBMC_datasets/raw'
adata_sets = []
for s in os.listdir(cr_dir):
    if not s.endswith('.csv'):
        sample_adata = sc.read_10x_h5(f'{cr_dir}/{s}/filtered_feature_bc_matrix.h5')
        #sample_adata = sc.read_10x_h5(cr_dir+"/"+s+"/"+"filtered_feature_bc_matrix.h5")
        sample_adata.var_names_make_unique()
        sample_adata.obs['library_id'] = s
        # scrub = scr.Scrublet(sample_adata.X)
        # doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
        # sample_adata.obs['doublet_scores'] = doublet_scores
        # sample_adata.obs['predicted_doublets'] = predicted_doublets
        adata_sets.append(sample_adata)
adata = adata_sets[0].concatenate(adata_sets[1:],join='outer')

adata.var['ensembl_id'] = adata.var['gene_ids']
adata.var['gene_ids'] = adata.var.index
adata.var['MT'] = adata.var_names.str.startswith('mt-')
adata.var['ribo'] = adata.var_names.str.startswith(('Rps','Rpl'))
sc.pp.calculate_qc_metrics(adata,qc_vars=['MT','ribo'],percent_top=[10,20],log1p=False,inplace=True)

adata.write_h5ad('10x_PBMC.h5ad')
