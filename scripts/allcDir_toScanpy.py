import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc
from fastparquet import write
path=sys.argv[1]
outpre=sys.argv[2]


files={'gmCH':[],'gmCG':[],'bmCH':[],'bmCG':[]}
for file in os.listdir(path):
    print(file)
    if 'gene_mCH' in file: files['gmCH'].append(file)
    elif 'gene_mCG' in file: files['gmCG'].append(file)
    elif 'bin_mCH' in file: files['bmCH'].append(file)
    elif 'bin_mCG' in file: files['bmCG'].append(file)

for ftype in files:
    mtype=''
    if 'CH' in ftype: mtype='H'
    else: mtype='G'
    #run first file of of loop
    df=pd.read_csv(path+files[ftype][0])
    df=df.rename(columns={"Unnamed: 0": "region"})
    df=df.set_index("region")
    df=df[~df['chrom'].isin(['chrY'])]

    obs=files[ftype][0][:-4]
    df[obs]=df['mC'+mtype]
    X_df=df[[obs]]

    covX_df=df[['C'+mtype]]
    covX_df=covX_df.rename(columns={'C'+mtype:obs})
    X_df
    count=0
    for file in files[ftype][1:]:
        if count%100==0: print(count)
        df=pd.read_csv(path+file)
        df=df.rename(columns={"Unnamed: 0": "region"})
        df=df.set_index("region")
        #df=df[~df['chrom'].isin(['chrX','chrY'])]
        df=df[~df['chrom'].isin(['chrY'])]

        obs=file[:-4]
        df[obs]=df['mC'+mtype]

        X_df[obs]=df[[obs]]
        covX_df[obs]=df[['C'+mtype]]


        count+=1
    adata=sc.AnnData(X_df.transpose())
    adata.write_h5ad(path+outpre+'_'+ftype+'.h5ad')
    print('Writing to: '+path+outpre+'_'+ftype+'xCov.parq')
    write(path+outpre+'_'+ftype+'xCov.parq',covX_df)

