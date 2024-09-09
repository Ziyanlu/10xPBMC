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

def get_markers(adata,groupby,key="rank_genes_groups",p_val_cutoff=0.05,logfc_cutoff=0.5):
    """\
    Extract markers from adata into Seurat-like table
    Extracts markers after they are computed by ``scanpy``. Produces Seurat-like
    table with fields
    ``"p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene"``
    Calculates the percentage of cells that express a given gene
    in the target cluster (``pct.1`` field) and outside the cluster
    (``pct.2`` field) from ``adata.raw`` matrix.
    Parameters
    ----------
    adata
        Annotated data matrix.
    groupby
        ``adata.obs`` field used for marker calculation
    key
        ``adata.uns`` key that has computed markers
    p_val_cutoff
        Drop all genes with adjusted p-value greater than or equal to this
    logfc_cutoff
        Drop all genes with average logFC less than or equal to this
    Returns
    -------
    Returns a pandas dataframe with above listed columns, optionally
    subsetted on the genes that pass the cutoffs.
    ``p_val`` field is a copy of adjusted p-value field.
    Example
    -------
    >>> sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", n_genes=200)
    >>> markers = sc_utils.get_markers(adata, "leiden")
    >>> markers.to_csv("markers.csv")
    """
    markers = pd.concat([
        pd.DataFrame(adata.uns[key]["names"]).melt(),
        pd.DataFrame(adata.uns[key]["pvals_adj"]).melt(),
        pd.DataFrame(adata.uns[key]["logfoldchanges"]).melt()
    ], axis=1)
    markers.columns = ("cluster", "gene", "cluster2", "p_val_adj", "cluster3", "avg_logFC")
    markers = markers.loc[:, ["cluster", "gene", "avg_logFC", "p_val_adj"]]
    markers = markers.loc[markers.avg_logFC > logfc_cutoff, ]
    markers = markers.loc[markers.p_val_adj < p_val_cutoff, ]
    markers["pct.1"] = pd.Series(dtype=float)
    markers["pct.2"] = pd.Series(dtype=float)
    for cluster in markers.cluster.unique():
        cells = adata.obs[groupby] == cluster
        in_cluster_selector = markers.cluster == cluster
        genes = markers.gene[in_cluster_selector]
        in_cluster = np.sum(adata.raw[cells, genes].X > 0, axis=0) / cells.sum()
        markers.loc[in_cluster_selector, "pct.1"] = in_cluster.T.A1
        other_cells = adata.obs[groupby] != cluster
        other_clusters = np.sum(adata.raw[other_cells, genes].X > 0, axis=0) / other_cells.sum()
        markers.loc[in_cluster_selector, "pct.2"] = other_clusters.T.A1
    markers["p_val"] = markers.p_val_adj
    markers = markers.loc[:, ["p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene"]]
    return markers
    
# cr_dir = '/projects/b1038/Pulmonary/sfenske/sequencing/data/2024-4/cellranger'
# adata_sets = []
# for s in os.listdir(cr_dir):
#     if not s.endswith('.csv'):
#         sample_adata = sc.read_10x_h5(
#             f'{cr_dir}/{s}/outs/filtered_feature_bc_matrix.h5')
#         sample_adata.var_names_make_unique()
#         sample_adata.obs['sample_id'] = s
#         sample_adata.var['MT'] = sample_adata.var_names.str.startswith('MT-')
#         sample_adata.var['ribo'] = sample_adata.var_names.str.startswith(('RPS','RPL'))
#         try:
#             scrub = scr.Scrublet(sample_adata.X)
#             doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=True)
#             sample_adata.obs['doublet_scores'] = doublet_scores
#             sample_adata.obs['predicted_doublets'] = predicted_doublets
#         except:
#             print(f"error for {s}",flush=True)
#             sample_adata.obs['doublet_scores'] = 0
#             sample_adata.obs['predicted_doublets'] = False
#         adata_sets.append(sample_adata)
# adata = adata_sets[0].concatenate(adata_sets[1:],join='outer')

# adata.var.index = [x.split('GRCh38______')[1] if 'GRCh38______' in x else x for x in adata.var.index]
# adata.var.gene_ids = [x.split('GRCh38______')[1] if 'GRCh38______' in x else x for x in adata.var.gene_ids]
# adata.var['ensembl_id'] = adata.var['gene_ids']
# adata.var['gene_ids'] = adata.var.index
# adata.var['MT'] = adata.var_names.str.startswith('MT-')
# adata.var['ribo'] = adata.var_names.str.startswith(('RPS','RPL'))
# sc.pp.calculate_qc_metrics(adata,qc_vars=['MT','ribo'],percent_top=[10,20],log1p=False,inplace=True)

# # for scArches
# adata.obs['sample'] = adata.obs['sample_id']

    
adata = sc.read_h5ad('../data/10x_PBMC.h5ad')

FOLDER = '../data'
OBJ_NAME = 'PBMC_integrated'
BATCH_VARIABLE = 'library_id'

adata.layers['counts'] = adata.X
sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata,subset=True,n_top_genes=2000,batch_key=BATCH_VARIABLE,flavor='seurat_v3',layer='counts')

param_set = {"n_layers":2,"n_hidden":256,"dropout_rate":0.2,"n_latent":10}
scvi.data.setup_anndata(adata,layer="counts",batch_key=BATCH_VARIABLE)
# vae = scvi.model.SCVI.load(f"../DeepSeq_1113_vae_2000", adata)

vae = scvi.model.SCVI(adata,**param_set)

train_kwargs = {
    "early_stopping": True,
    "early_stopping_monitor": "reconstruction_loss_validation",
    "early_stopping_patience": 25
}
vae.train(max_epochs=400,**train_kwargs)
vae.save(f"../code/{OBJ_NAME}",overwrite=True)


latent_adata = vae.get_latent_representation()
adata.obsm["X_scVI"] = latent_adata
adata.layers["scvi_normalized"] = vae.get_normalized_expression(library_size=10e4)
sc.pp.neighbors(adata,use_rep='X_scVI')
sc.tl.umap(adata)
sc.tl.leiden(adata,key_added="leiden_scVI")

adata.uns['log1p']["base"] = None
sc.tl.rank_genes_groups(adata, "leiden_scVI", method="wilcoxon", n_genes=200)
markers = get_markers(adata,'leiden_scVI')

for col in adata.obs.columns:
    if adata.obs[col].dtype == 'category':
        adata.obs[col] = adata.obs[col].astype('str')


if not os.path.isdir(FOLDER):
    os.mkdir(FOLDER)
markers.to_csv(f'{FOLDER}/{OBJ_NAME}-markers.csv')
adata.obs.to_csv(f'{FOLDER}/{OBJ_NAME}-metadata.csv')
adata.write_h5ad(f'{FOLDER}/{OBJ_NAME}.h5ad')
