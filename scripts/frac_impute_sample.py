import os 
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import fastparquet as fp

sample_root=sys.argv[1]
min_mean_feature_threshold=sys.argv[2]
feature_impute_threshold=sys.argv[3]

class methyl_matrix:
    def __init__(self,counts_file,coverage_file):
        self.adata=sc.read_h5ad(counts_file)
        self.cov=fp.ParquetFile(coverage_file).to_pandas()
        self.cov=self.cov.transpose()
        self.cov=self.cov.rename(columns={i:str(i) for i in self.cov.columns})
    def load_metadata(self,metadata_file):
        df=pd.read_csv(metadata_file,delimiter='\t').set_index('Sample')
        for col in df.columns:
            self.adata.obs[col]=df[col]
    def filter_features(self,coverage_threshold):
        print(len(self.adata.var.index))
        mean_cov=self.cov.describe().transpose()['mean']
        self.adata.var['mean_coverage']=mean_cov
        use_features=self.adata.var[self.adata.var['mean_coverage']>coverage_threshold].index
        print(len(use_features))
        self.adata=self.adata[:,use_features]
        self.cov=self.cov[use_features]
        
    def calc_mc_frac(self):
        self.adata.X=self.adata.X/self.cov
        
    def impute(self,min_cov):
        self.adata.var['mean_X']=[np.ma.array(i, mask=np.isnan(i)).mean() for i in self.adata.X.T]
        self.adata.var['mean_X']=self.adata.var['mean_X'].astype(float)
        pass_qc=self.cov>min_cov
        X=self.adata.to_df()
        for col in X:
            meanX=self.adata.var.loc[col]['mean_X']
            X[col]=np.where(pass_qc[col],X[col],meanX)
        self.adata.X=X

sr=sample_root
for mod in ['gmCG','gmCH','bmCG','bmCH']:
    print(mod)
    mc_mat=methyl_matrix('sc_files/'+sr+'_'+mod+'.h5ad','sc_files/'+sr+'_'+mod+'xCov.parq')
    mc_mat.adata.to_df().to_csv('mc_mats/'+sr+'_'+mod+'_counts.tsv.gz',sep='\t',compression='gzip')
    mc_mat.cov.to_csv('mc_mats/'+sr+'_'+mod+'_coverage.tsv.gz',sep='\t',compression='gzip')
    mc_mat.filter_features(min_mean_feature_threshold)
    mc_mat.calc_mc_frac()
    mc_mat.impute(feature_impute_threshold)
    mc_mat.load_metadata('metadata_files/'+sr+'.summary.txt')
    mc_mat.adata.write_h5ad('mc_mats/'+sr+'_'+mod+'.h5ad')
