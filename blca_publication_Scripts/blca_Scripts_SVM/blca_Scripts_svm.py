import anndata
import pandas as pd

import cstarpy
import os
import numpy as np
from cstarpy.separation import CellStateTransition

datadir = '/home/jing/Phd_project/project_UCD_blca/blca_publication_OUTPUT/'
os.chdir(datadir)

mtx =os.path.join('blca_publication_OUTPUT_sct/',"sct_corrected_UMI.mtx")
cells=pd.read_csv(os.path.join('blca_publication_OUTPUT_sct/',"sct_corrected_UMI_cells.txt"),header=None,index_col=0)
features=pd.read_csv(os.path.join('blca_publication_OUTPUT_sct/','sct_corrected_UMI_genes.txt'),header=None,sep='\t',index_col=0)
adata = sc.read_mtx(mtx)

adata.obs['CellID']= cells.index.tolist()
adata.var['Gene']= features.index.tolist()
adata.var.index= adata.var['Gene']
display(adata)
print(f'Max val of sct matrix before transformation, {np.max(adata.X)}')

m_h_conversion = pd.read_csv('/home/jing/Phd_project/project_UCD_blca/blca_OUTPUT/blca_OUTPUT_m_h_convert/m_h_convertion.csv',index_col='Mouse_Genes')
common_mouse_genes= list(set(adata.var.index).intersection(set(m_h_conversion.index)))

# Iterate over common genes and use .at for proper assignment
adata.var['Lincs'] = np.nan
for i in common_mouse_genes:
    if i in m_h_conversion.index:
        # Use .at for efficient and correct value assignment
        adata.var.at[i, 'Lincs'] = m_h_conversion.at[i, 'Human_Genes']

adata.obs['Sample'] = adata.obs['CellID'].str.split('_').str[0]
adata.obs.index = adata.obs['CellID']

adata.obs['Type'] = 'N/A'
for i in ['GSM5288668', 'GSM5288669']:
    adata.obs.loc[adata.obs['Sample'].str.contains(i, na=False), 'Type'] = 'NMIBC'
for i in ['GSM5288670', 'GSM5288671']:
    adata.obs.loc[adata.obs['Sample'].str.contains(i, na=False), 'Type'] = 'MIBC'    
for i in ['GSM5288672', 'GSM5288674']:
    adata.obs.loc[adata.obs['Sample'].str.contains(i, na=False), 'Type'] = 'Healthy'    

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata,base=2)
print(f'Max val after log1p, {np.max(adata.to_df())}')

adata_lincs = adata[:, adata.var['Lincs'].notnull()]
MIBC= adata_lincs[adata_lincs.obs['Type'] == 'MIBC']

mibc_df= MIBC.to_df()


NMIBC= adata_lincs[adata_lincs.obs['Type'] == 'NMIBC']
nmibc_df = NMIBC.to_df()

#
healthy= adata_lincs[adata_lincs.obs['Type'] == 'Healthy']
healthy_df = healthy.to_df()


merged=adata_lincs[adata_lincs.obs['Type']!='N/A']
merged=merged[merged.obs['Type']!='Healthy']
merged_df = merged.to_df()


cst_onc = CellStateTransition('onc',healthy_df, merged_df)
dpd_onc_scores=cst_onc.get_dpd()
dpd_onc_scores

norm_s_onc_df = pd.DataFrame(np.stack([cst_onc.n, cst_onc.s], axis=1), 
                             index=cst_onc.svm_input.data.columns, columns=["n", "s"])

#Run
cst = CellStateTransition('test', nmibc_df, mibc_df)
dpd_scores = cst.get_dpd()

norm_s_df = pd.DataFrame(np.stack([cst.n, cst.s], axis=1), index=cst.svm_input.data.columns, columns=["n", "s"])

pd.to_pickle(dpd_scores,"dpd_lincs_sct_inv_blca_N_M.pkl")
pd.to_pickle(norm_s_df,"/stv_lincs_sct_inv_blca_N_M.pkl")

pd.to_pickle(dpd_onc_scores,"dpd_lincs_sct_onc_healthy_onc.pkl")
pd.to_pickle(norm_s_onc_df,"stv_lincs_sct_onc_healthy_onc.pkl")

data_all_cells = cst_onc.svm_input.data
DPD_invasive = cst.get_dpd(data_all_cells)
DPD_onc = cst_onc.get_dpd(data_all_cells)
DPD_all = pd.concat([DPD_invasive, DPD_onc], axis=1)

pd.to_pickle(DPD_all,"plot_lincs_sct_onc_inv.pkl")
