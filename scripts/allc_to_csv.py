import os
import numpy as np
import pandas as pd
import gzip
import sys

allc_file_path=sys.argv[1]
bin_size=int(sys.argv[2])
reference_genome_path=sys.argv[3]
genome_annotation_bed=sys.argv[4]
out_prefix=sys.argv[5]


def import_genome_annotations(bedfile):
    adf=pd.read_csv(bedfile,compression='gzip',sep='\t',header=None)
    adf=adf.rename(columns={0:'chrom',1:'start',2:'end',3:'name'})
    adf=adf.set_index('name')
    return adf


def make_bin_df(bin_size,chr_lens):
    bin_regions_dict={}
    bin_count=0
    for key in chr_lens.keys():
        chr_bins=[(i,i+bin_size-1) for i in range(0,chr_lens[key],bin_size)]
        chr_bins[-1]=(chr_bins[-1][0],chr_lens[key]-1)
        for b in chr_bins:
            bin_regions_dict[bin_count]={'chrom':key.split(' ')[0],'start':b[0],'end':b[1]}
            bin_count+=1
    bdf=pd.DataFrame(bin_regions_dict).transpose()
    bdf.index.name='name'
    return bdf


def get_region_mC(allc_df,start,stop):
    allc_df=allc_df[allc_df[1].between(start, stop, inclusive=False)]
    cg_df=allc_df[allc_df[3].isin(['CGA','CGT','CGG','CGC'])]
    ch_df=allc_df[~allc_df[3].isin(['CGA','CGT','CGG','CGC'])]
    cg_df=cg_df[cg_df[5]<=2]
    ch_df=ch_df[ch_df[5]<=2]
    return({"mCH":ch_df[4].sum(),"CH":ch_df[5].sum(),"mCG":cg_df[4].sum(),"CG":cg_df[5].sum()})
   
   
#this could be more efficient
#reduce runtime by using pandas interface rather than for i,r in iterrows in nthis function
def process_allc_file(filename,region_df):
    allc_df=pd.read_csv(filename,compression='gzip',sep='\t',header=None)
    allc_df=allc_df[~allc_df[0].isin(['chrY'])]
    ch_region_df={}
    cg_region_df={}
    for chrm in allc_df[0].unique():
        chrm_allc_df=allc_df[allc_df[0]==chrm]
        cr_df=region_df[region_df['chrom']==chrm]
        for i,r in cr_df.iterrows():
            region_info=get_region_mC(chrm_allc_df,r['start'],r['end'])
            entry={"chrom":chrm,"start":r['start'],"end":r['end']}
            entry["mCH"]=region_info["mCH"]
            entry["CH"]=region_info["CH"]
            ch_region_df[i]=entry.copy()
            entry["mCG"]=region_info["mCG"]
            entry["CG"]=region_info["CG"]
            cg_region_df[i]=entry.copy()
    return {'ch_df':pd.DataFrame(ch_region_df),'cg_df':pd.DataFrame(cg_region_df)}

chr_len_dict={}
with open(reference_genome_path,'r') as f:
    line=f.readline()[:-1]
    while line:
        spline=line.split('\t')
        chr_len_dict[spline[0]]=int(spline[1])
        line=f.readline()[:-1]

region_df=import_genome_annotations(genome_annotation_bed)
region_df['start']=region_df['start']-2000
#region_df['start']=region_df['start']
region_df['end']=region_df['end']+2000
#region_df['end']=region_df['end']
region_df=region_df[~region_df['chrom'].isin(['chrY'])]

bin_df=make_bin_df(bin_size,chr_len_dict)
#bin_df=bin_df[~bin_df['chrom'].isin(['chrX','chrY'])]

outname=out_prefix+allc_file_path.split('/')[-1].split('.')[0]
print(outname)
rm_dfs=process_allc_file(allc_file_path,region_df)
rm_dfs['ch_df'].transpose().to_csv(outname+'_gene_mCH.csv')
rm_dfs['cg_df'].transpose().to_csv(outname+'_gene_mCG.csv')
bm_dfs=process_allc_file(allc_file_path,bin_df)
bm_dfs['ch_df'].transpose().to_csv(outname+'_bin_mCH.csv')
bm_dfs['cg_df'].transpose().to_csv(outname+'_bin_mCG.csv')
